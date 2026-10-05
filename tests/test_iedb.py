"""Tests for the IEDB epitope-id lookup.

These exercise the candidate-selection rule offline, against the cached responses, so the
suite makes no network calls. The rule matters: ranking by assay count instead of specificity
gave E03 and E06 the same id, because storage-protein sequences are repetitive enough that one
registered peptide can contain several panel cores.
"""

from __future__ import annotations

import csv

import pytest

from grains import config, iedb


def test_best_candidate_prefers_the_shortest_containing_peptide():
    rows = [
        {"structure_id": 462152, "linear_sequence": "PQPEQPFPQPEQPFPQPEQPFPQP", "n_assays": 9},
        {"structure_id": 2265304, "linear_sequence": "QPFPQPEQPFPQS", "n_assays": 1},
    ]
    # The 24-mer has far more assays; the 13-mer is the more specific identifier.
    assert iedb.best_candidate(rows)["structure_id"] == 2265304


def test_best_candidate_breaks_length_ties_by_assay_count_then_id():
    rows = [
        {"structure_id": 300, "linear_sequence": "PFPQPELPYAA", "n_assays": 1},
        {"structure_id": 100, "linear_sequence": "PFPQPELPYBB", "n_assays": 5},
        {"structure_id": 200, "linear_sequence": "PFPQPELPYCC", "n_assays": 5},
    ]
    assert iedb.best_candidate(rows)["structure_id"] == 100


def test_best_candidate_handles_no_rows():
    assert iedb.best_candidate([]) is None


@pytest.mark.skipif(not iedb.OUT_PATH.exists(), reason="IEDB lookup not run yet")
def test_shared_ids_are_flagged_not_hidden():
    """An id reused across rows must carry a warning; it identifies neither row alone."""
    rows = list(csv.DictReader(iedb.OUT_PATH.open(newline="")))
    by_id: dict[str, list[dict]] = {}
    for row in rows:
        if row["iedb_id"]:
            by_id.setdefault(row["iedb_id"], []).append(row)
    for iedb_id, sharing in by_id.items():
        if len(sharing) > 1:
            for row in sharing:
                assert row["note"].startswith("WARNING"), (
                    f"{row['epitope_id']} shares IEDB id {iedb_id} but carries no warning"
                )


@pytest.mark.skipif(not iedb.OUT_PATH.exists(), reason="IEDB lookup not run yet")
def test_resolved_ids_reference_a_peptide_containing_the_core():
    """Whatever entry a row points at must actually contain that row's core."""
    for row in csv.DictReader(iedb.OUT_PATH.open(newline="")):
        if row["iedb_id"] and row["core_9mer"]:
            assert row["core_9mer"] in row["iedb_sequence"], row["epitope_id"]


@pytest.mark.skipif(not iedb.OUT_PATH.exists(), reason="IEDB lookup not run yet")
def test_peptides_csv_iedb_column_matches_the_lookup():
    recorded = {
        row["epitope_id"]: row["iedb_id"]
        for row in csv.DictReader(config.PEPTIDES_CSV.open(newline=""))
    }
    for row in csv.DictReader(iedb.OUT_PATH.open(newline="")):
        assert recorded[row["epitope_id"]] == row["iedb_id"], row["epitope_id"]


def test_e03_and_e06_cores_genuinely_overlap():
    """Documents why those two rows share an id: 7 of 9 residues in common.

    If this ever stops holding, the shared id would be a lookup bug rather than biology,
    and the warning in peptides_NOTES.md would need revisiting.
    """
    rows = {row["epitope_id"]: row for row in
            csv.DictReader(config.PEPTIDES_CSV.open(newline=""))}
    e03, e06 = rows["E03"]["core_9mer"], rows["E06"]["core_9mer"]
    overlaps = [e03[i:] for i in range(1, 9) if e03[i:] == e06[: 9 - i]]
    assert overlaps == ["PQPEQPF"]
    assert len(overlaps[0]) == 7
