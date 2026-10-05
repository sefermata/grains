"""Tests for the peptide-MHC superposition / RMSD baseline.

The regression these guard is the 4OZI chain-pairing trap: that structure's peptide chain J
binds MHC A/B while chain I binds C/D, the reverse of the chain lettering. Pairing by
lettering yields a clean-looking groove fit and a ~62 A peptide RMSD -- a wrong number that
arrives as a result rather than an error.
"""

from __future__ import annotations

import pytest

from grains import superpose
from grains.superpose import ComplexSpec

pytestmark = pytest.mark.skipif(
    not (superpose.PDB_DIR / "1S9V.pdb").exists()
    or not (superpose.PDB_DIR / "4OZI.pdb").exists(),
    reason="validation PDB files not downloaded (data/raw/pdb is git-ignored)",
)


def test_matched_pairs_identical_sequences():
    pairs = superpose.matched_pairs("PFPQPELPY", "PFPQPELPY")
    assert pairs == [(i, i) for i in range(9)]


def test_matched_pairs_skips_mismatched_columns():
    # Native vs deamidated E01 core differ at exactly one position (Q->E at P6).
    pairs = superpose.matched_pairs("PFPQPQLPY", "PFPQPELPY")
    assert (5, 5) not in pairs
    assert len(pairs) == 8


def test_matched_pairs_handles_offset_and_tag():
    # 1S9V's peptide carries an extra N-terminal L; 4OZI's carries a C-terminal PGS tag.
    pairs = superpose.matched_pairs("LQPFPQPELPY", "QPFPQPELPYPGS")
    assert len(pairs) == 10
    assert pairs[0] == (1, 0)


def test_core_window_found_and_missing():
    assert superpose.core_window("LQPFPQPELPY", superpose.E01_CORE) == list(range(2, 11))
    assert superpose.core_window("QQYPSGEGSFQPSQENPQ", superpose.E01_CORE) is None


def test_wrong_peptide_chain_is_rejected():
    """4OZI A/B + I is the plausible-looking but wrong pairing; it must raise, not score."""
    wrong = ComplexSpec("4OZI_wrong", superpose.PDB_DIR / "4OZI.pdb", "A", "B", "I")
    with pytest.raises(ValueError, match="not bound to"):
        superpose.compare(wrong, superpose.S9V_COPY1)


def test_correct_pairings_pass_verification():
    from Bio.PDB import PDBParser

    parser = PDBParser(QUIET=True)
    for spec in (superpose.S9V_COPY1, superpose.S9V_COPY2,
                 superpose.OZI_COPY1, superpose.OZI_COPY2):
        structure = parser.get_structure(spec.label, spec.pdb_path)
        distances = superpose.verify_pairing(structure, spec)
        # A bound peptide sits within a few Angstroms of both chains of the groove.
        assert max(distances.values()) < 10.0


def test_baseline_rmsds_are_physically_plausible():
    """A good groove fit forces a small peptide RMSD -- the peptide is inside the groove."""
    rows = superpose.run()
    assert len(rows) == 4
    for row in rows:
        assert row["mhc_fit_rmsd_A"] < 2.0
        assert row["peptide_rmsd_A"] < 3.0, row
        assert row["core_rmsd_A"] is not None
        assert row["n_core_ca"] == 9


def test_same_crystal_copies_are_tighter_than_cross_deposition():
    """The noise floor must be tighter than the cross-deposition ceiling, or the
    baseline tells us nothing about what a model should be held to."""
    rows = {row["comparison"]: row for row in superpose.run()}
    same_crystal = rows["1S9V_ABC -> 1S9V_DEF"]["core_rmsd_A"]
    cross = rows["1S9V_ABC -> 4OZI_ABJ"]["core_rmsd_A"]
    assert same_crystal < cross
