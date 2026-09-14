"""Build the teaching dataset for the Data 8 grains project.

Reads the seed-storage-protein FASTA files fetched for the main project
(../../data/raw/*.fasta) and writes three CSVs into ./data/:

    grains.csv              labeled training/exploration set
    grains_mystery.csv      unlabeled hold-out set (students classify these)
    grains_mystery_key.csv  answers for grains_mystery.csv (instructor only)

Design choices (kept simple on purpose):
  * only the 20 standard amino acids; sequences with X/B/Z/U or '*' are dropped
  * sequence length restricted to 50-600 residues (drops fragments and outliers)
  * exact-duplicate sequences removed within a grain
  * a fixed RNG seed (8, naturally) makes the sample reproducible
  * NO composition columns are written - computing them is Part 1 of the project

Run:  python build_dataset.py
"""

from __future__ import annotations

import csv
import random
import re
from pathlib import Path

RAW = Path(__file__).resolve().parents[2] / "data" / "raw"
OUT = Path(__file__).resolve().parent / "data"

STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")
MIN_LEN, MAX_LEN = 50, 600
SEED = 8

# fasta basename stem -> (grain label, gluten yes/no)
SOURCES = {
    "bread-wheat__gliadin-glutenin": ("wheat", "yes"),
    "barley__hordein": ("barley", "yes"),
    "rye__secalin": ("rye", "yes"),
    "rice__glutelin-orzenin": ("rice", "no"),
    "maize-corn__zein": ("maize", "no"),
    "sorghum__kafirin": ("sorghum", "no"),
    "oat__avenin": ("oat", "no"),
    "quinoa__11s-2s": ("quinoa", "no"),
    "foxtail-millet__setarin": ("foxtail_millet", "no"),
    "buckwheat__13s": ("buckwheat", "no"),
}

# how many sequences to keep per grain in the labeled set (capped at availability)
LABELED_CAP = {
    "wheat": 80, "barley": 999, "rye": 999,
    "rice": 25, "maize": 25, "sorghum": 20, "oat": 999,
    "quinoa": 20, "foxtail_millet": 20, "buckwheat": 20,
}

# how many sequences to pull aside (before sampling) into the mystery set
MYSTERY_TAKE = {"wheat": 5, "barley": 4, "rice": 4, "maize": 3, "oat": 2}


def parse_fasta(path: Path):
    """Yield (accession, protein_name, sequence) from one of our FASTA files."""
    acc = name = None
    seq: list[str] = []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if acc is not None:
                yield acc, name, "".join(seq)
            head = line[1:]
            before_os = head.split(" OS=")[0]
            id_block, _, rest = before_os.partition(" ")
            acc = id_block.split("|")[0]
            name = rest.strip() or "unknown"
            seq = []
        else:
            seq.append(line.strip())
    if acc is not None:
        yield acc, name, "".join(seq)


def clean_records(path: Path) -> list[dict]:
    seen: set[str] = set()
    out: list[dict] = []
    for acc, name, seq in parse_fasta(path):
        seq = seq.upper()
        if not (MIN_LEN <= len(seq) <= MAX_LEN):
            continue
        if set(seq) - STANDARD_AA:
            continue
        if seq in seen:
            continue
        seen.add(seq)
        out.append(
            {
                "accession": acc,
                "protein_name": re.sub(r"\s+", " ", name)[:80],
                "sequence": seq,
            }
        )
    return out


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    rng = random.Random(SEED)

    labeled: list[dict] = []
    mystery: list[dict] = []

    for stem, (grain, gluten) in SOURCES.items():
        path = RAW / f"{stem}.fasta"
        if not path.exists():
            print(f"  skip (missing): {path.name}")
            continue
        recs = clean_records(path)
        rng.shuffle(recs)

        n_mystery = MYSTERY_TAKE.get(grain, 0)
        for r in recs[:n_mystery]:
            mystery.append({**r, "grain": grain, "gluten": gluten})

        pool = recs[n_mystery:]
        cap = LABELED_CAP.get(grain, 999)
        for r in pool[:cap]:
            labeled.append({**r, "grain": grain, "gluten": gluten})

        print(f"  {grain:15} available={len(recs):4}  mystery={n_mystery}  labeled={min(cap, len(pool))}")

    rng.shuffle(labeled)
    rng.shuffle(mystery)

    for i, r in enumerate(labeled, 1):
        r["seq_id"] = f"S{i:03d}"
    for i, r in enumerate(mystery, 1):
        r["seq_id"] = f"M{i:02d}"

    cols = ["seq_id", "grain", "gluten", "protein_name", "accession", "sequence"]
    _write(OUT / "grains.csv", labeled, cols)
    _write(OUT / "grains_mystery.csv", mystery, ["seq_id", "protein_name", "accession", "sequence"])
    _write(OUT / "grains_mystery_key.csv", mystery, ["seq_id", "grain", "gluten"])

    n_yes = sum(r["gluten"] == "yes" for r in labeled)
    print(
        f"\ngrains.csv: {len(labeled)} rows "
        f"({n_yes} gluten / {len(labeled) - n_yes} gluten-free) "
        f"| baseline accuracy = {max(n_yes, len(labeled) - n_yes) / len(labeled):.1%}"
    )
    print(f"grains_mystery.csv: {len(mystery)} rows")


def _write(path: Path, rows: list[dict], cols: list[str]) -> None:
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {path.relative_to(OUT.parent)}  ({len(rows)} rows)")


if __name__ == "__main__":
    main()
