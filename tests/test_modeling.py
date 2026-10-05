"""Tests for the two modeling routes and the scoring restriction between them.

The important test here is the licensing guard: an AlphaFold Server model must never reach
a binding-affinity scorer, and a PRODIGY number computed from one would look identical to
a valid number. Keeping that as a test means the constraint survives refactoring.
"""

from __future__ import annotations

import json

import pytest

from grains import modeling


def test_only_pandora_may_be_scored():
    modeling.assert_scorable(modeling.ROUTE_PANDORA)  # must not raise


def test_alphafold_route_is_refused_for_scoring():
    with pytest.raises(ValueError, match="Prohibited Use Policy"):
        modeling.assert_scorable(modeling.ROUTE_ALPHAFOLD, context="PRODIGY dG")


def test_unknown_route_is_refused_for_scoring():
    """Default-deny: a route nobody has vetted is not scorable just because it's new."""
    with pytest.raises(ValueError):
        modeling.assert_scorable("some_new_folding_service")


def test_validation_peptide_matches_the_crystal_structure_peptide():
    """The model and the crystal must be the same molecule or the RMSD needs a caveat."""
    from grains import superpose

    assert modeling.VALIDATION_PEPTIDE == "LQPFPQPELPY"
    assert superpose.E01_CORE in modeling.VALIDATION_PEPTIDE


def test_deamidated_core_carries_a_glutamate():
    """TG2 deamidation is Q->E; E is standard, so no PTM is declared to AlphaFold Server.
    If this peptide ever lost its E the submission would silently model the native form."""
    assert "E" in modeling.VALIDATION_PEPTIDE
    core_start = modeling.VALIDATION_PEPTIDE.find("PFPQPELPY")
    assert modeling.VALIDATION_PEPTIDE[core_start + 5] == "E"


@pytest.mark.skipif(
    not modeling.HLA_FASTA.exists(),
    reason="hla_dq25.fasta not generated yet (run grains.fetch_hla)",
)
def test_submission_json_has_the_expected_shape(tmp_path):
    out = modeling.write_validation_submission(tmp_path / "job.json")
    jobs = json.loads(out.read_text())

    assert isinstance(jobs, list) and len(jobs) == 1
    job = jobs[0]
    assert job["dialect"] == "alphafoldserver"
    assert job["version"] == 1
    assert job["modelSeeds"] == []

    # Three chains: MHC alpha, MHC beta, peptide.
    sequences = [entry["proteinChain"]["sequence"] for entry in job["sequences"]]
    assert len(sequences) == 3
    assert all(entry["proteinChain"]["count"] == 1 for entry in job["sequences"])
    assert sequences[2] == modeling.VALIDATION_PEPTIDE

    alpha, beta = modeling.load_hla_chains()
    assert sequences[0] == alpha
    assert sequences[1] == beta


@pytest.mark.skipif(
    not modeling.HLA_FASTA.exists(),
    reason="hla_dq25.fasta not generated yet (run grains.fetch_hla)",
)
def test_hla_chains_are_the_1s9v_constructs():
    alpha, beta = modeling.load_hla_chains()
    assert (len(alpha), len(beta)) == (193, 198)
    assert "X" not in alpha and "X" not in beta


def test_pandora_availability_reports_a_reason_when_unavailable():
    ok, detail = modeling.pandora_available()
    assert isinstance(ok, bool)
    assert detail
    if not ok:
        # An unavailability message has to say what to do about it.
        assert "modeller" in detail.lower() or "pandora" in detail.lower()
