"""Extract the HLA-DQ2.5 alpha/beta ectodomain sequences for modeling.

Rather than pulling full-length chains from UniProt and guessing where to trim the
signal peptide / transmembrane / cytoplasmic segments, this reads them straight out of
the chosen Checkpoint 2 validation structure (data/raw/pdb/1S9V.pdb, see
data/reference/pdb_structures.md). Those chains are already exactly the groove+Ig
ectodomain construct the crystallographers solved for DQA1*05:01 (chain A) and
DQB1*02:01 (chain B) -- no domain-boundary guessing, and it guarantees the sequences
used for AlphaFold Server / PANDORA modeling match the structure they'll be superposed
against for the RMSD check.

Reads the `SEQRES` records (the full designed construct), not the resolved `ATOM`
coordinates: chain B has a genuine 7-residue gap (a disordered loop the crystallographers
couldn't resolve, PDB numbering 106-112) where the two sides of the gap are NOT adjacent
in the actual structure. Walking CA atoms and concatenating whatever's present would
silently splice across that gap into a sequence that skips 7 real residues with no sign
anything is missing -- wrong input for AlphaFold Server, which needs the true designed
sequence, disordered residues included.

Usage
-----
    PYTHONPATH=src python -m grains.fetch_hla

Output
------
    data/processed/hla_dq25.fasta   two records: DQA1*05:01 and DQB1*02:01 ectodomains
"""

from __future__ import annotations

from Bio import SeqIO

from . import config

SOURCE_PDB = config.RAW / "pdb" / "1S9V.pdb"
ALPHA_CHAIN = "A"
BETA_CHAIN = "B"


def chain_sequences(pdb_path) -> dict[str, str]:
    """Chain id -> full designed sequence, from SEQRES (not the resolved ATOM records)."""
    seqs: dict[str, str] = {}
    for record in SeqIO.parse(pdb_path, "pdb-seqres"):
        chain_id = record.id.split(":")[-1]
        seqs[chain_id] = str(record.seq)
    return seqs


def run() -> dict[str, str]:
    config.ensure_dirs()
    if not SOURCE_PDB.exists():
        raise FileNotFoundError(
            f"{SOURCE_PDB} not found -- download it first "
            "(curl https://files.rcsb.org/download/1S9V.pdb -o data/raw/pdb/1S9V.pdb)"
        )
    chains = chain_sequences(SOURCE_PDB)
    alpha_seq = chains[ALPHA_CHAIN]
    beta_seq = chains[BETA_CHAIN]

    records = {
        "DQA1*05:01_ectodomain|source=1S9V_chainA": alpha_seq,
        "DQB1*02:01_ectodomain|source=1S9V_chainB": beta_seq,
    }

    out_path = config.PROCESSED / "hla_dq25.fasta"
    lines: list[str] = []
    for header, seq in records.items():
        lines.append(f">{header} len={len(seq)}")
        lines.extend(seq[i : i + 60] for i in range(0, len(seq), 60))
    out_path.write_text("\n".join(lines) + "\n")

    print(f"wrote {out_path}")
    for header, seq in records.items():
        n_x = seq.count("X")
        flag = f"  ({n_x} non-standard residue(s) — check before modeling)" if n_x else ""
        print(f"  {header}: {len(seq)} aa{flag}")

    return records


if __name__ == "__main__":
    run()
