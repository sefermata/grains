"""Tests for the PANDORA route: the Biopython/IMGT compatibility layer and anchor derivation.

All offline. Nothing here builds a model or touches the network -- the point is the two places
the PANDORA route can go silently wrong:

1. The compatibility shims. Two of the three are stubs for code paths this project does not
   use, and they must *raise* rather than return something plausible if a future change
   wanders onto them.
2. Anchor derivation. PANDORA does not validate allele names or anchor positions; it scores
   every template and returns the best. Wrong anchors therefore produce a model, not an error.
"""

from __future__ import annotations

import csv
import shutil

import pytest

from grains import config, modeling, pandora_compat


# --- compatibility shims ----------------------------------------------------------------

def test_pandora_imports_with_the_shims():
    """The whole point of pandora_compat: PANDORA 0.9 is importable under Biopython 1.88."""
    pytest.importorskip("PANDORA")
    from PANDORA.Database import Database  # noqa: F401
    from PANDORA.Pandora import Pandora  # noqa: F401
    from PANDORA.PMHC import PMHC  # noqa: F401


def test_removed_biopython_names_are_present_again():
    from Bio.PDB import Polypeptide

    assert hasattr(Polypeptide, "three_to_one")
    from Bio import pairwise2

    assert hasattr(pairwise2.align, "globalxx")
    from Bio.Align.Applications import MuscleCommandline  # noqa: F401


def test_three_to_one_stub_raises_rather_than_guessing():
    """Off this project's path. If something reaches it, that must be loud, not approximate."""
    from Bio.PDB import Polypeptide

    with pytest.raises(NotImplementedError, match="protein_letters_3to1"):
        Polypeptide.three_to_one("ALA")


def test_pairwise2_stub_raises_rather_than_guessing():
    from Bio import pairwise2

    with pytest.raises(NotImplementedError, match="PairwiseAligner"):
        pairwise2.align.globalxx("AAA", "AAA")


def test_muscle_shim_emits_v3_syntax():
    """PANDORA does `os.system(str(msl) + ' -quiet')`, so __str__ is the entire contract."""
    rendered = str(pandora_compat.MuscleShim(input="/tmp/in.fasta", out="/tmp/out.afa"))
    assert " -in " in rendered and " -out " in rendered
    assert "/tmp/in.fasta" in rendered and "/tmp/out.afa" in rendered
    # v5's flags would mean MUSCLE exits 0 and writes nothing.
    assert "-align" not in rendered and "-output" not in rendered


def test_muscle_shim_quotes_paths_with_spaces():
    rendered = str(pandora_compat.MuscleShim(input="/tmp/a b/in.fasta", out="/tmp/out.afa"))
    assert "'/tmp/a b/in.fasta'" in rendered


def test_muscle_shim_explains_itself_when_muscle_is_absent(monkeypatch):
    monkeypatch.setattr(pandora_compat.shutil, "which", lambda _: None)
    with pytest.raises(FileNotFoundError, match="muscle=3.8.1551"):
        str(pandora_compat.MuscleShim(input="in", out="out"))


@pytest.mark.skipif(shutil.which("muscle") is None, reason="MUSCLE not installed")
def test_installed_muscle_speaks_the_syntax_the_shim_emits():
    ok, detail = pandora_compat.muscle_syntax_ok()
    assert ok, detail


def test_imgt_query_url_is_https():
    """PANDORA hardcodes http, and imgt.org serves nothing over plain HTTP."""
    assert pandora_compat.IMGT_QUERY_URL.startswith("https://")


# --- anchor derivation ------------------------------------------------------------------

def test_anchors_for_the_validation_peptide():
    """E01 in 1S9V: core PFPQPELPY starts at position 3 of LQPFPQPELPY."""
    anchors = modeling.class_ii_anchors(modeling.VALIDATION_PEPTIDE, "PFPQPELPY")
    assert anchors == [3, 6, 8, 11]


def test_anchor_p6_is_the_deamidated_glutamate_for_e01():
    """Ties the geometry back to the biology: TG2's Q->E lands on a charge-preferring pocket."""
    peptide = modeling.VALIDATION_PEPTIDE
    anchors = modeling.class_ii_anchors(peptide, "PFPQPELPY")
    p6 = anchors[2]
    assert peptide[p6 - 1] == "E"


def test_anchors_rejects_a_core_that_is_not_nine_residues():
    with pytest.raises(ValueError, match="9 residues"):
        modeling.class_ii_anchors("PQPEQPFPWQPEQ", "PQPEQPF")


def test_anchors_rejects_a_core_absent_from_the_extended_sequence():
    """The E09 failure mode from Checkpoint 1: a deamidated core against a native window."""
    with pytest.raises(ValueError, match="not a substring"):
        modeling.class_ii_anchors("QGSFQPSQQNPQ", "EGSFQPSQE")


def _cored_dq25_rows() -> list[dict]:
    with config.PEPTIDES_CSV.open(newline="") as handle:
        return [row for row in csv.DictReader(handle)
                if len(row["core_9mer"]) == 9 and row["restriction"] == "DQ2.5"]


def test_every_panel_row_yields_in_range_anchors():
    """Anchors are derived per row, never hardcoded -- cores sit at different offsets."""
    rows = _cored_dq25_rows()
    assert rows, "no cored DQ2.5 rows found"
    for row in rows:
        peptide = row["extended_sequence"]
        anchors = modeling.class_ii_anchors(peptide, row["core_9mer"])
        assert len(anchors) == 4, row["epitope_id"]
        assert anchors == sorted(set(anchors)), row["epitope_id"]
        assert all(1 <= position <= len(peptide) for position in anchors), row["epitope_id"]


def test_cores_sit_at_different_offsets_so_a_fixed_anchor_list_would_be_wrong():
    """Guards the reason class_ii_anchors exists at all."""
    offsets = {row["extended_sequence"].find(row["core_9mer"]) for row in _cored_dq25_rows()}
    assert len(offsets) > 1, f"all cores at the same offset {offsets}; derivation is moot"


def test_deamidated_dq25_cores_carry_a_glutamate_at_a_charge_pocket():
    """DQ2.5's P4/P6/P7 pockets are positively charged and want a negative side chain.

    Note P7 is included: it is a charge-preferring pocket for DQ2.5 but is *not* one of the
    four geometric anchors PANDORA uses (P1/P4/P6/P9). E08 (PQQSFPEQQ) is the panel row where
    that distinction bites -- its only glutamate is at P7.
    """
    for row in _cored_dq25_rows():
        if row["form"] != "deamidated":
            continue
        core = row["core_9mer"]
        positions = [index + 1 for index, residue in enumerate(core) if residue == "E"]
        assert any(position in (4, 6, 7) for position in positions), (
            f"{row['epitope_id']} ({core}) has glutamate(s) at {positions}, none at a "
            "DQ2.5 charge-preferring pocket (P4/P6/P7)"
        )


def test_dq25_allele_names_match_pandoras_imgt_spelling():
    """PANDORA stores `HLA-DQA1*0501`: asterisk kept, four digits, no colon.

    This test exists because the first guess was `HLA-DQA10501` and it was wrong. The
    spelling is not derivable from the standard nomenclature, so it is pinned here.
    """
    assert modeling.DQ25_ALLELES == ["HLA-DQA1*0501", "HLA-DQB1*0201"]
    for name in modeling.DQ25_ALLELES:
        assert "*" in name and ":" not in name
        assert name.startswith("HLA-DQ")


# --- template selection (needs the built database) --------------------------------------

def _database():
    from grains import pandora_db

    ready, _ = pandora_db.database_ready()
    if not ready:
        pytest.skip("PANDORA template database not built")
    return modeling.load_database()


def test_dq25_allele_names_are_real_names_in_the_database():
    modeling.check_allele_names(modeling.DQ25_ALLELES, _database())


def test_a_misspelled_allele_is_rejected_with_near_misses():
    database = _database()
    with pytest.raises(ValueError, match="not present"):
        modeling.check_allele_names(["HLA-DQA10501"], database)   # the wrong first guess


def test_dq25_alleles_admit_the_known_dq25_templates():
    """The allele filter must leave a real candidate set, not one structure."""
    eligible = modeling.templates_for_alleles(modeling.DQ25_ALLELES, _database())
    assert "1S9V" in eligible
    assert len(eligible) >= 10, eligible


def test_dqb1_0201_alone_would_collapse_to_1s9v():
    """The trap this spelling choice avoids.

    `HLA-DQB1*0201` is a correct, specific name for DQ2.5's beta chain, but PANDORA parsed it
    from only one structure's IMGT remarks. Requesting it alone is therefore a hard filter
    down to 1S9V, and *every* panel peptide would be modeled on the Checkpoint 2 validation
    structure -- silently, since PANDORA reports no problem. Pairing it with HLA-DQA1*0501 is
    what keeps the candidate set real.
    """
    eligible = modeling.templates_for_alleles(["HLA-DQB1*0201"], _database())
    assert eligible == ["1S9V"], eligible


def test_the_validation_structure_is_a_template():
    """Checkpoint 3's gate exists because of this. If it ever stops being true, re-read
    validate_pandora's three conditions -- `no_exact` would become meaningless."""
    assert "1S9V" in _database().MHCII_data


def test_gate_conditions_each_retain_eligible_templates():
    from grains import validate_pandora

    database = _database()
    for name, spec in validate_pandora.conditions(database).items():
        subset = validate_pandora.restricted(database, spec["withhold"])
        eligible = modeling.templates_for_alleles(modeling.DQ25_ALLELES, subset)
        assert eligible, f"condition {name} leaves PANDORA no template"


def test_no_core_withholds_every_template_sharing_the_e01_core():
    from grains import validate_pandora
    from grains.superpose import E01_CORE

    database = _database()
    spec = validate_pandora.conditions(database)["no_core"]
    subset = validate_pandora.restricted(database, spec["withhold"])
    for pdb_id in modeling.templates_for_alleles(modeling.DQ25_ALLELES, subset):
        assert E01_CORE not in (subset.MHCII_data[pdb_id].peptide or ""), pdb_id


def test_longest_shared_run_is_the_longest_contiguous_match():
    from grains.validate_pandora import longest_shared_run

    assert longest_shared_run("LQPFPQPELPY", "QPFPQPELPYP") == "QPFPQPELPY"
    assert longest_shared_run("AAAA", "BBBB") == ""
    assert longest_shared_run("ABCDE", "XCDEY") == "CDE"


# --- template conformational state ------------------------------------------------------

def test_tcr_detection_separates_known_binary_from_known_ternary():
    """1S9V is a bare peptide-MHC complex; 4OZI has a T-cell receptor bound.

    Checkpoint 1 established that all of 4OZF-4OZI are TCR-bound, which is why 1S9V and not
    4OZF became the validation structure. This re-derives it from IMGT's own annotation.
    """
    from grains import validate_pandora

    if validate_pandora.is_tcr_bound("1S9V") is None:
        pytest.skip("IMGT flat files not on disk")
    assert validate_pandora.is_tcr_bound("1S9V") is False
    assert validate_pandora.is_tcr_bound("4OZI") is True
    assert validate_pandora.is_tcr_bound("6MFG") is False


def test_most_dq25_templates_are_tcr_bound():
    """The reason template state has to be controlled rather than assumed.

    8 of the 12 DQ2.5 templates are ternary, so left alone PANDORA will often pick one and
    inherit a peptide conformation from a complex that includes a TCR -- while this project's
    v1 scope is deliberately peptide-MHC without a TCR.
    """
    from grains import validate_pandora

    database = _database()
    if validate_pandora.is_tcr_bound("1S9V") is None:
        pytest.skip("IMGT flat files not on disk")
    dq25 = modeling.templates_for_alleles(modeling.DQ25_ALLELES, database)
    ternary = validate_pandora.tcr_bound(database)
    binary = sorted(set(dq25) - ternary)
    assert binary == ["1S9V", "6MFF", "6MFG", "6U3M"], binary
    assert len(set(dq25) & ternary) == 8


def test_binary_condition_leaves_a_similar_template_available():
    """The controlled comparison only works if a binary template of comparable similarity
    survives -- 6MFG shares 10 of 11 residues with the target, exactly as 4OZI does."""
    from grains import validate_pandora

    database = _database()
    if validate_pandora.is_tcr_bound("1S9V") is None:
        pytest.skip("IMGT flat files not on disk")
    spec = validate_pandora.conditions(database)["no_exact_binary"]
    subset = validate_pandora.restricted(database, spec["withhold"])
    eligible = modeling.templates_for_alleles(modeling.DQ25_ALLELES, subset)
    assert "6MFG" in eligible
    shared = validate_pandora.longest_shared_run(
        modeling.VALIDATION_PEPTIDE, subset.MHCII_data["6MFG"].peptide)
    assert len(shared) == 10, shared


def test_binary_only_database_excludes_every_ternary_template():
    from grains import pandora_db

    database = _database()
    if pandora_db.is_tcr_bound("1S9V") is None:
        pytest.skip("IMGT flat files not on disk")
    subset = modeling.binary_templates_only(database)
    assert not (set(subset.MHCII_data) & pandora_db.tcr_bound(database))
    assert len(subset.MHCII_data) < len(database.MHCII_data)
    # The original must not be mutated -- validate_pandora compares against the full set.
    assert "4OZI" in database.MHCII_data


def test_binary_only_still_leaves_dq25_templates_to_model_on():
    database = _database()
    from grains import pandora_db

    if pandora_db.is_tcr_bound("1S9V") is None:
        pytest.skip("IMGT flat files not on disk")
    eligible = modeling.templates_for_alleles(
        modeling.DQ25_ALLELES, modeling.binary_templates_only(database))
    assert eligible == ["1S9V", "6MFF", "6MFG", "6U3M"], eligible
