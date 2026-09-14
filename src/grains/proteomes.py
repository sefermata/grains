"""Fetch seed-storage proteins for each grain from the UniProt REST API.

Usage
-----
    python -m grains.proteomes                 # all grains in data/reference/grains.csv
    python -m grains.proteomes --grain wheat   # one grain (substring match on common_name)
    python -m grains.proteomes --force         # ignore cached FASTA and re-download

Outputs
-------
    data/raw/<grain>__<family>.fasta           # one FASTA per grain
    data/processed/proteome_index.csv          # accession-level metadata for every entry
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import pandas as pd
import requests

from . import config


def _slug(text: str) -> str:
    keep = [c if c.isalnum() else "-" for c in text.lower()]
    return "".join(keep).strip("-").replace("--", "-")


def _family_key(family_field: str) -> str:
    """Map a grains.csv storage_protein_family value to a STORAGE_PROTEIN_TERMS key."""
    token = family_field.split()[0].lower().strip()
    if token in config.STORAGE_PROTEIN_TERMS:
        return token
    # fall back to a generic prolamin / seed-storage search
    return "prolamin"


def build_query(taxon_id: int, terms: list[str]) -> str:
    name_clause = " OR ".join(f'protein_name:"{t}"' for t in terms)
    return f"(taxonomy_id:{taxon_id}) AND ({name_clause})"


def fetch_entries(query: str, *, page_size: int | None = None) -> list[dict]:
    """Return a list of simplified UniProt entry dicts for a query (handles paging)."""
    page_size = page_size or config.UNIPROT_PAGE_SIZE
    params = {
        "query": query,
        "fields": "accession,id,protein_name,organism_name,length,reviewed,sequence",
        "format": "json",
        "size": str(page_size),
    }
    headers = {"User-Agent": f"grains-project ({config.CONTACT_EMAIL})"}
    url = config.UNIPROT_BASE
    out: list[dict] = []

    while url:
        resp = requests.get(
            url, params=params, headers=headers, timeout=config.REQUEST_TIMEOUT
        )
        resp.raise_for_status()
        payload = resp.json()
        for entry in payload.get("results", []):
            out.append(_simplify(entry))
        # UniProt puts the next page in the Link header; params only on first call
        url = resp.links.get("next", {}).get("url")
        params = None
        if url:
            time.sleep(0.2)  # be polite
    return out


def _simplify(entry: dict) -> dict:
    names = entry.get("proteinDescription", {})
    rec = names.get("recommendedName", {}).get("fullName", {}).get("value")
    sub = None
    if not rec:
        subs = names.get("submissionNames", [])
        if subs:
            sub = subs[0].get("fullName", {}).get("value")
    return {
        "accession": entry.get("primaryAccession"),
        "entry_name": entry.get("uniProtkbId"),
        "protein_name": rec or sub or "",
        "organism": entry.get("organism", {}).get("scientificName", ""),
        "length": entry.get("sequence", {}).get("length"),
        "reviewed": entry.get("entryType", "").lower().startswith("uniprotkb reviewed"),
        "sequence": entry.get("sequence", {}).get("value", ""),
    }


def _write_fasta(records: list[dict], path: Path) -> int:
    lines: list[str] = []
    n = 0
    for r in records:
        seq = r["sequence"]
        if not seq:
            continue
        if not (config.MIN_SEQ_LEN <= len(seq) <= config.MAX_SEQ_LEN):
            continue
        header = (
            f">{r['accession']}|{r['entry_name']} {r['protein_name']} "
            f"OS={r['organism']} len={len(seq)} "
            f"reviewed={'yes' if r['reviewed'] else 'no'}"
        )
        lines.append(header)
        lines.extend(seq[i : i + 60] for i in range(0, len(seq), 60))
        n += 1
    path.write_text("\n".join(lines) + ("\n" if lines else ""))
    return n


def _parse_fasta_headers(path: Path) -> list[dict]:
    """Recover per-entry rows from a FASTA we wrote earlier (see _write_fasta header format)."""
    rows: list[dict] = []
    for line in path.read_text().splitlines():
        if not line.startswith(">"):
            continue
        head = line[1:]
        acc_block, _, rest = head.partition(" ")
        accession, _, entry_name = acc_block.partition("|")
        length = reviewed = None
        for tok in rest.split():
            if tok.startswith("len="):
                length = int(tok[4:])
            elif tok.startswith("reviewed="):
                reviewed = tok[9:] == "yes"
        protein_name = rest.split(" OS=")[0].strip()
        rows.append(
            {
                "accession": accession,
                "entry_name": entry_name,
                "protein_name": protein_name,
                "length": length,
                "reviewed": reviewed,
            }
        )
    return rows


def run(grain_filter: str | None = None, *, force: bool = False) -> pd.DataFrame:
    config.ensure_dirs()
    grains = pd.read_csv(config.GRAINS_CSV)
    if grain_filter:
        mask = grains["common_name"].str.contains(grain_filter, case=False, na=False)
        grains = grains[mask]
        if grains.empty:
            sys.exit(f"no grain matches {grain_filter!r} in {config.GRAINS_CSV}")

    entry_rows: list[dict] = []
    summary_rows: list[dict] = []

    for _, g in grains.iterrows():
        fam_key = _family_key(str(g["storage_protein_family"]))
        terms = config.STORAGE_PROTEIN_TERMS[fam_key]
        out_path = config.RAW / f"{_slug(g['common_name'])}__{_slug(fam_key)}.fasta"

        if out_path.exists() and not force:
            entries = _parse_fasta_headers(out_path)
            source = "cached"
            print(f"  cached  {out_path.name:<40} {len(entries)} seqs")
        else:
            query = build_query(int(g["ncbi_taxon_id"]), terms)
            print(f"  fetch   {g['common_name']:<16} taxon={g['ncbi_taxon_id']} terms={terms}")
            try:
                hits = fetch_entries(query)
            except requests.HTTPError as exc:
                print(f"    ! HTTP error for {g['common_name']}: {exc}")
                continue
            kept = _write_fasta(hits, out_path)
            print(f"    -> {len(hits)} hits, {kept} kept after length filter -> {out_path.name}")
            entries = [
                {
                    "accession": e["accession"],
                    "entry_name": e["entry_name"],
                    "protein_name": e["protein_name"],
                    "length": len(e["sequence"]),
                    "reviewed": e["reviewed"],
                }
                for e in hits
                if e["sequence"]
                and config.MIN_SEQ_LEN <= len(e["sequence"]) <= config.MAX_SEQ_LEN
            ]
            source = "fetched"

        for e in entries:
            entry_rows.append(
                {
                    "common_name": g["common_name"],
                    "celiac_status": g["celiac_status"],
                    "family": fam_key,
                    **e,
                }
            )
        summary_rows.append(
            {
                "common_name": g["common_name"],
                "species": g["species"],
                "celiac_status": g["celiac_status"],
                "family": fam_key,
                "search_terms": "|".join(terms),
                "n_entries": len(entries),
                "n_reviewed": sum(1 for e in entries if e["reviewed"]),
                "fasta": out_path.name if entries else "",
                "source": source,
            }
        )

    index = pd.DataFrame(entry_rows)
    summary = pd.DataFrame(summary_rows)
    if not index.empty:
        idx_path = config.PROCESSED / "proteome_index.csv"
        index.to_csv(idx_path, index=False)
        print(f"\nwrote {idx_path} ({len(index)} protein entries)")
    if not summary.empty:
        sum_path = config.PROCESSED / "proteome_summary.csv"
        summary.to_csv(sum_path, index=False)
        print(f"wrote {sum_path} ({len(summary)} grains)")
        low = summary[summary["n_entries"] < 5]["common_name"].tolist()
        if low:
            print(f"\n  LOW COVERAGE (<5 entries): {', '.join(low)} — note in the report")
    return index


def _parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--grain", help="substring of common_name to fetch just one grain")
    p.add_argument("--force", action="store_true", help="re-download even if cached")
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    args = _parse_args(argv)
    run(args.grain, force=args.force)


if __name__ == "__main__":
    main()
