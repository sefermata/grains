"""Attach an IEDB epitope ID to every panel peptide (Checkpoint 1, checklist step 2).

The IEDB (Immune Epitope Database) is the community registry of experimentally observed
epitopes. Recording its ID per row makes each peptide traceable to the assay evidence behind
it, independently of the Sollid nomenclature the panel was built from -- a second, orthogonal
identifier for the same molecule.

Two lookups, in order, because neither alone covers the panel:

1. **Exact match on `core_9mer`.** Some epitopes are registered as the bare 9-mer core
   (`PFPQPELPY` is IEDB 47539). When this hits it is the cleanest possible answer.
2. **Contains-search on the core.** Most are registered as the longer peptide that was
   actually synthesised and assayed, with the core inside it -- E09's core sits in eight
   registered peptides, among them `QQYPSGEGSFQPSQENPQ`, which is exactly 2NNA's deposited
   peptide. Exact search on our `extended_sequence` finds nothing for any row, since our
   window is chosen from the gene rather than copied from whatever length a given lab used.

Where the contains-search returns several candidates, the one with the most assay records is
chosen (most experimental support), ties broken by lowest IEDB id so the choice is
deterministic. **Every candidate id is kept** in `iedb_candidate_ids` -- picking one for the
reference table shouldn't silently discard the others.

Responses are cached on disk, so re-runs cost no requests and the column can be rebuilt
offline. Delete the cache file to force a refresh.

Usage
-----
    PYTHONPATH=src python -m grains.iedb

Output
------
    data/processed/peptide_iedb_map.csv   full lookup result, all candidates
    data/reference/peptides.csv           `iedb_id` column added/updated
"""

from __future__ import annotations

import csv
import json
import time
from pathlib import Path

import requests

from . import config

QUERY_API = "https://query-api.iedb.org/epitope_search"
CACHE_PATH = config.PROCESSED / ".iedb_cache.json"
OUT_PATH = config.PROCESSED / "peptide_iedb_map.csv"
REQUEST_PAUSE = 0.6     # seconds between calls; the API is free and public, so be gentle
CONTAINS_TIMEOUT = 180  # a contains-search scans every registered sequence; it needs room
RETRIES = 3

FIELDNAMES = [
    "epitope_id", "canonical_name", "core_9mer", "iedb_id", "iedb_sequence",
    "iedb_match_type", "iedb_n_assays", "iedb_n_candidates", "iedb_candidate_ids", "note",
]


def _load_cache() -> dict:
    if CACHE_PATH.exists():
        return json.loads(CACHE_PATH.read_text())
    return {}


def _save_cache(cache: dict) -> None:
    CACHE_PATH.write_text(json.dumps(cache, indent=0, sort_keys=True))


def query(sequence: str, mode: str, cache: dict | None = None) -> list[dict]:
    """Search IEDB for `sequence`. `mode` is "exact" or "contains".

    Returns the raw result rows, trimmed to the fields used here. Cached by (mode, sequence).
    """
    cache = cache if cache is not None else {}
    key = f"{mode}:{sequence}"
    if key in cache:
        return cache[key]

    operator = f"eq.{sequence}" if mode == "exact" else f"like.*{sequence}*"
    # A contains-search is a scan over every registered epitope sequence and routinely runs
    # well past the 30 s used for UniProt, so it gets its own budget and a couple of retries.
    # An exact search is indexed and returns promptly.
    timeout = config.REQUEST_TIMEOUT if mode == "exact" else CONTAINS_TIMEOUT
    last_error: Exception | None = None
    for attempt in range(RETRIES):
        try:
            response = requests.get(
                QUERY_API,
                params={"linear_sequence": operator, "limit": 50},
                headers={"User-Agent": f"grains-project ({config.CONTACT_EMAIL})"},
                timeout=timeout,
            )
            response.raise_for_status()
            break
        except (requests.Timeout, requests.ConnectionError) as exc:
            last_error = exc
            if attempt == RETRIES - 1:
                raise
            time.sleep(2 ** attempt)
    rows = [
        {
            "structure_id": row.get("structure_id"),
            "linear_sequence": row.get("linear_sequence"),
            "n_assays": len(row.get("iedb_assay_ids") or []),
        }
        for row in response.json()
    ]
    cache[key] = rows
    time.sleep(REQUEST_PAUSE)
    return rows


def best_candidate(rows: list[dict]) -> dict | None:
    """Shortest registered peptide containing the core; then most assays; then lowest id.

    Shortest-first, not best-studied-first. The id is meant to identify *this* epitope, and
    storage-protein sequences are highly repetitive: IEDB 462152 is a 24-residue omega-gliadin
    repeat (`PQPEQPFPQPEQPFPQPEQPFPQP`) that contains the cores of both E03 and E06. Ranking
    by assay count picked that one entry for two different epitopes -- an id shared between
    rows identifies neither. The shortest containing peptide is the most specific match, so
    it is the one that earns the row. Assay count then breaks ties toward the better-evidenced
    entry, and the id makes the choice deterministic.
    """
    if not rows:
        return None
    return sorted(
        rows,
        key=lambda r: (len(r["linear_sequence"]), -r["n_assays"], r["structure_id"]),
    )[0]


def lookup(core: str, cache: dict) -> dict:
    """Resolve one core to an IEDB entry, preferring an exact 9-mer registration."""
    if not core:
        return {"iedb_match_type": "no_core", "candidates": []}

    exact = query(core, "exact", cache)
    if exact:
        return {"iedb_match_type": "exact_core", "candidates": exact,
                "chosen": best_candidate(exact)}

    contains = query(core, "contains", cache)
    if contains:
        return {"iedb_match_type": "contains_core", "candidates": contains,
                "chosen": best_candidate(contains)}

    return {"iedb_match_type": "not_found", "candidates": []}


def run(write_reference: bool = True) -> list[dict]:
    config.ensure_dirs()
    cache = _load_cache()

    with config.PEPTIDES_CSV.open(newline="") as handle:
        peptides = list(csv.DictReader(handle))
        fields = list(peptides[0].keys())

    rows: list[dict] = []
    try:
        for peptide in peptides:
            core = (peptide.get("core_9mer") or "").strip()
            result = lookup(core, cache)
            chosen = result.get("chosen")
            candidates = result["candidates"]

            note = ""
            if result["iedb_match_type"] == "no_core":
                note = "row has no core_9mer (33-mer context row); not looked up"
            elif result["iedb_match_type"] == "not_found":
                note = "core not registered in IEDB, alone or inside a longer peptide"
            elif result["iedb_match_type"] == "contains_core" and len(candidates) > 1:
                note = (f"core occurs in {len(candidates)} registered peptides; chose the one "
                        "with the most assay records (all ids kept in iedb_candidate_ids)")

            rows.append({
                "epitope_id": peptide["epitope_id"],
                "canonical_name": peptide["canonical_name"],
                "core_9mer": core,
                "iedb_id": chosen["structure_id"] if chosen else "",
                "iedb_sequence": chosen["linear_sequence"] if chosen else "",
                "iedb_match_type": result["iedb_match_type"],
                "iedb_n_assays": chosen["n_assays"] if chosen else "",
                "iedb_n_candidates": len(candidates),
                "iedb_candidate_ids": ";".join(str(c["structure_id"]) for c in candidates),
                "note": note,
            })
    finally:
        _save_cache(cache)

    # An id reused across two panel rows identifies neither of them. The repetitive
    # storage-protein sequences make this a live hazard, so flag it in the output rather
    # than leaving two rows quietly pointing at the same entry.
    seen: dict[str, list[str]] = {}
    for row in rows:
        if row["iedb_id"]:
            seen.setdefault(str(row["iedb_id"]), []).append(row["epitope_id"])
    shared = {iedb_id: ids for iedb_id, ids in seen.items() if len(ids) > 1}
    for row in rows:
        others = shared.get(str(row["iedb_id"]), [])
        if others:
            row["note"] = (
                f"WARNING: IEDB id {row['iedb_id']} is shared with {', '.join(o for o in others if o != row['epitope_id'])}"
                " -- the registered peptide contains more than one panel core, so it does not "
                "identify this epitope uniquely. " + row["note"]
            ).strip()

    with OUT_PATH.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDNAMES)
        writer.writeheader()
        writer.writerows(rows)

    if write_reference:
        if "iedb_id" not in fields:
            fields.insert(fields.index("verified"), "iedb_id")
        by_id = {row["epitope_id"]: row for row in rows}
        for peptide in peptides:
            peptide["iedb_id"] = by_id[peptide["epitope_id"]]["iedb_id"]
        with config.PEPTIDES_CSV.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields)
            writer.writeheader()
            writer.writerows(peptides)

    markers = {"exact_core": "exact", "contains_core": "within",
               "not_found": "NOT FOUND", "no_core": "-"}
    for row in rows:
        marker = markers[row["iedb_match_type"]]
        n_candidates = row["iedb_n_candidates"]
        extra = f"({n_candidates} candidates)" if n_candidates > 1 else ""
        print(f"  {row['epitope_id']:>4} {row['canonical_name']:<26} "
              f"{str(row['iedb_id']):>8}  {marker:<9} "
              f"{row['iedb_sequence']:<20} {extra}")

    found = sum(1 for r in rows if r["iedb_id"])
    cored = sum(1 for r in rows if r["core_9mer"])
    print(f"\nwrote {OUT_PATH}")
    print(f"{found}/{cored} cored rows resolved to an IEDB epitope id")
    return rows


if __name__ == "__main__":
    run()
