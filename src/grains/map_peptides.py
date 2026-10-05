"""Locate each panel peptide in a real UniProt sequence: accession + residue range.

Checkpoint 1, checklist step 3 ("confirm the extended_sequence is a real substring of a
UniProt entry for that protein + grain, record accession + residue range") -- with a wrinkle
that step 3 as written doesn't survive contact with the data:

**A deamidated peptide is not a substring of any native proteome.** That's the whole point of
deamidation: TG2 converts specific glutamines to glutamate, so `QLQPFPQPELPYPQPQL` (E01) has
an E where the gene encodes a Q. Searching for it literally finds nothing, and concluding "not
in UniProt" would be wrong.

So the search allows a peptide `E` to match a proteome `E` *or* `Q`, and records which
positions had to flex. That turns one check into two: it both locates the peptide (step 3) and
independently confirms the deamidation sites (step 4) -- a position that had to flex from Q is
a deamidation site, and it should be one the literature already claims. Other residues must
match exactly, so this is not a fuzzy search: for E01 there is exactly one degree of freedom
per E, and a 17-residue match with everything else literal is not landing by chance.

Reported per peptide:
  * `match_type` -- `exact` (found literally, no flexing: a native-form peptide, or a
    deamidated one whose E is genuinely encoded) or `deamidated` (needed >=1 Q->E).
  * `deamidated_positions` -- 1-based offsets within the peptide where proteome Q became E.
  * `found_in_declared_grain` -- whether the hit came from the grain the row claims. A peptide
    found only in a different grain is a data error worth seeing, not something to paper over.

Misses are reported, not hidden. The rye/oat/teff proteome coverage gap is real and known
(`data/processed/proteome_summary.csv`): rye has 3 sequences and oat 2, so rye/oat epitopes may
genuinely have no local match. That's a gap in public annotation, not a bad epitope.

Usage
-----
    PYTHONPATH=src python -m grains.map_peptides

Output
------
    data/processed/peptide_uniprot_map.csv
"""

from __future__ import annotations

import csv
import re
from dataclasses import dataclass
from pathlib import Path

from Bio import SeqIO

from . import config

# peptides.csv's `source_grain` is colloquial ("wheat"); the proteome fastas are per
# grains.csv common_name ("bread wheat", "durum wheat"). Wheat epitopes may sit in either
# the hexaploid or the tetraploid entries, so both files are searched.
GRAIN_FASTAS: dict[str, list[str]] = {
    "wheat": ["bread-wheat__gliadin-glutenin.fasta", "durum-wheat__gliadin-glutenin.fasta"],
    "barley": ["barley__hordein.fasta"],
    "rye": ["rye__secalin.fasta"],
    "oat": ["oat__avenin.fasta"],
}

OUT_PATH = config.PROCESSED / "peptide_uniprot_map.csv"

FIELDNAMES = [
    "epitope_id", "canonical_name", "source_grain", "source_protein", "form",
    "extended_sequence", "uniprot_accession", "uniprot_entry_name", "protein_name",
    "start", "end", "match_type", "n_deamidation_sites", "deamidated_positions",
    "found_in_declared_grain", "matched_grain", "n_other_hits", "note",
]


@dataclass(frozen=True)
class Protein:
    accession: str
    entry_name: str
    description: str
    sequence: str
    grain: str


def load_proteins(grain_files: dict[str, list[str]] | None = None) -> dict[str, list[Protein]]:
    """Grain -> proteins, parsed from the fastas written by `grains.proteomes`."""
    grain_files = grain_files or GRAIN_FASTAS
    proteins: dict[str, list[Protein]] = {}
    for grain, filenames in grain_files.items():
        entries: list[Protein] = []
        for filename in filenames:
            path = config.RAW / filename
            if not path.exists():
                continue
            for record in SeqIO.parse(path, "fasta"):
                # Header: ">P10388|GLT5_WHEAT Glutenin, ... OS=... len=848 reviewed=yes"
                accession, _, rest = record.id.partition("|")
                entry_name = rest or ""
                description = record.description.split(" OS=")[0]
                description = description.partition(" ")[2] or description
                entries.append(
                    Protein(accession, entry_name, description, str(record.seq), grain)
                )
        proteins[grain] = entries
    return proteins


def deamidation_tolerant_pattern(peptide: str) -> re.Pattern[str]:
    """Compile `peptide` with every E allowed to match E or Q.

    One degree of freedom per glutamate, every other residue literal. Q is NOT allowed to
    match E: deamidation goes Q->E in the peptide relative to the gene, never the reverse.
    """
    return re.compile("".join(f"[EQ]" if residue == "E" else re.escape(residue)
                              for residue in peptide))


def find_in_protein(peptide: str, protein: Protein) -> tuple[int, list[int]] | None:
    """(1-based start, deamidated offsets) for the first hit, or None.

    Prefers a literal match: if the peptide occurs verbatim the gene already encodes those
    glutamates and nothing was deamidated, which is a stronger result than a flexed match.
    """
    literal = protein.sequence.find(peptide)
    if literal >= 0:
        return literal + 1, []

    match = deamidation_tolerant_pattern(peptide).search(protein.sequence)
    if match is None:
        return None
    start = match.start()
    flexed = [
        offset + 1
        for offset, residue in enumerate(peptide)
        if residue == "E" and protein.sequence[start + offset] == "Q"
    ]
    return start + 1, flexed


def locate(peptide: str, declared_grain: str,
           proteins: dict[str, list[Protein]]) -> list[tuple[Protein, int, list[int]]]:
    """All hits for `peptide`, declared grain first so the best-provenance hit wins."""
    ordered = [declared_grain] + [g for g in proteins if g != declared_grain]
    hits: list[tuple[Protein, int, list[int]]] = []
    for grain in ordered:
        for protein in proteins.get(grain, []):
            found = find_in_protein(peptide, protein)
            if found is not None:
                hits.append((protein, found[0], found[1]))
    return hits


def run() -> list[dict]:
    config.ensure_dirs()
    proteins = load_proteins()
    counts = {grain: len(entries) for grain, entries in proteins.items()}
    print("proteins loaded per grain:", counts)
    print()

    rows: list[dict] = []
    with config.PEPTIDES_CSV.open(newline="") as handle:
        peptides = list(csv.DictReader(handle))

    for peptide_row in peptides:
        peptide = peptide_row["extended_sequence"]
        grain = peptide_row["source_grain"]
        hits = locate(peptide, grain, proteins)

        row: dict[str, object] = {
            "epitope_id": peptide_row["epitope_id"],
            "canonical_name": peptide_row["canonical_name"],
            "source_grain": grain,
            "source_protein": peptide_row["source_protein"],
            "form": peptide_row["form"],
            "extended_sequence": peptide,
        }

        if not hits:
            # Before blaming proteome coverage, check whether the *core* is present. If it
            # is, the peptide's own verified core sits in a real gene and only the flanking
            # residues of extended_sequence don't -- a curation error in the flanks, not a
            # missing sequence. Reporting "coverage gap" here would blame the data for a
            # problem in the panel, and send the next person looking in the wrong place.
            core = (peptide_row.get("core_9mer") or "").strip()
            core_hits = locate(core, grain, proteins) if core else []
            coverage = counts.get(grain, 0)

            if core_hits:
                protein, start, flexed = core_hits[0]
                row.update({
                    "uniprot_accession": protein.accession,
                    "uniprot_entry_name": protein.entry_name,
                    "protein_name": protein.description,
                    "start": start,
                    "end": start + len(core) - 1,
                    "match_type": "core_only",
                    "n_deamidation_sites": len(flexed),
                    "deamidated_positions": ";".join(str(p) for p in flexed),
                    "found_in_declared_grain": "yes" if protein.grain == grain else "no",
                    "matched_grain": protein.grain,
                    "n_other_hits": len(core_hits) - 1,
                    "note": (
                        f"extended_sequence not found in any gene, but core {core} IS present "
                        f"in {protein.accession} ({protein.entry_name}) at {start}-"
                        f"{start + len(core) - 1}. The flanking residues of extended_sequence "
                        "are therefore not gene-faithful -- a curation error in the flanks, "
                        "NOT a proteome coverage gap. extended_sequence is what gets modeled, "
                        "so this needs fixing before modeling."
                    ),
                })
            else:
                row.update({
                    "uniprot_accession": "", "uniprot_entry_name": "", "protein_name": "",
                    "start": "", "end": "", "match_type": "not_found",
                    "n_deamidation_sites": "", "deamidated_positions": "",
                    "found_in_declared_grain": "no", "matched_grain": "", "n_other_hits": 0,
                    "note": (
                        f"neither extended_sequence nor core found; {grain} has only "
                        f"{coverage} sequence(s) locally (known public-annotation "
                        "coverage gap)"
                        if coverage <= 5 else
                        "neither extended_sequence nor core found in any loaded proteome"
                    ),
                })
            rows.append(row)
            continue

        protein, start, flexed = hits[0]
        row.update({
            "uniprot_accession": protein.accession,
            "uniprot_entry_name": protein.entry_name,
            "protein_name": protein.description,
            "start": start,
            "end": start + len(peptide) - 1,
            "match_type": "exact" if not flexed else "deamidated",
            "n_deamidation_sites": len(flexed),
            "deamidated_positions": ";".join(str(p) for p in flexed),
            "found_in_declared_grain": "yes" if protein.grain == grain else "no",
            "matched_grain": protein.grain,
            "n_other_hits": len(hits) - 1,
            "note": "" if protein.grain == grain else (
                f"WARNING: not found in declared grain {grain!r}; "
                f"matched {protein.grain!r} instead -- check the row"
            ),
        })
        rows.append(row)

    with OUT_PATH.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDNAMES)
        writer.writeheader()
        writer.writerows(rows)

    for row in rows:
        status = row["match_type"]
        if status in ("not_found", "core_only"):
            label = "NOT FOUND" if status == "not_found" else "CORE ONLY"
            print(f"  {row['epitope_id']:>4} {row['canonical_name']:<26} {label}  "
                  f"({row['note']})")
        else:
            flag = "" if row["found_in_declared_grain"] == "yes" else "  <-- WRONG GRAIN"
            sites = row["deamidated_positions"] or "-"
            print(f"  {row['epitope_id']:>4} {row['canonical_name']:<26} "
                  f"{row['uniprot_accession']:<8} {row['start']:>4}-{row['end']:<4} "
                  f"{status:<11} Q->E at {sites}{flag}")

    tally: dict[str, int] = {}
    for row in rows:
        tally[str(row["match_type"])] = tally.get(str(row["match_type"]), 0) + 1

    print(f"\nwrote {OUT_PATH}")
    full = tally.get("exact", 0) + tally.get("deamidated", 0)
    print(f"{full}/{len(rows)} peptides fully located (whole extended_sequence in a gene): "
          f"{tally.get('deamidated', 0)} needed Q->E, {tally.get('exact', 0)} matched literally")
    if tally.get("core_only"):
        print(f"{tally['core_only']} rows are CORE ONLY -- verified core sits in a real gene "
              "but the extended_sequence flanks do not. These must be fixed before modeling; "
              "extended_sequence is the modeling input.")
    if tally.get("not_found"):
        print(f"{tally['not_found']} rows not found at all (proteome coverage gap)")
    return rows


if __name__ == "__main__":
    run()
