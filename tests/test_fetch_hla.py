"""Regression test for fetch_hla.py's SEQRES-based extraction.

1S9V chain B has a real 7-residue gap unresolved in the crystal (PDB numbering
106-112). Extracting by walking CA atoms instead of SEQRES silently splices across
that gap -- this test pins the full SEQRES-derived lengths so that regression can't
come back unnoticed.
"""

from __future__ import annotations

import pytest

from grains import config
from grains.fetch_hla import SOURCE_PDB, chain_sequences

pytestmark = pytest.mark.skipif(
    not SOURCE_PDB.exists(),
    reason=f"{SOURCE_PDB} not downloaded (data/raw/ is git-ignored; see checkpoints.md)",
)


def test_chain_lengths_match_full_seqres_not_resolved_atoms_only():
    chains = chain_sequences(SOURCE_PDB)
    # full designed construct lengths (from SEQRES), not the resolved-atom counts
    # (179 / 181) a CA-walk would silently produce by splicing across chain B's gap
    assert len(chains["A"]) == 193
    assert len(chains["B"]) == 198


def test_chain_sequences_are_plausible_protein():
    chains = chain_sequences(SOURCE_PDB)
    aa = set("ACDEFGHIKLMNPQRSTVWYX")
    for chain_id, seq in chains.items():
        assert set(seq) <= aa, f"bad residues in chain {chain_id}: {seq!r}"
        assert seq.count("X") == 0, f"unexpected non-standard residues in chain {chain_id}"
