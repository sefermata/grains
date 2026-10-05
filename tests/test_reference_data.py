"""Smoke tests for the reference data and package wiring."""

from __future__ import annotations

import pandas as pd
import pytest

from grains import config
from grains.proteomes import build_query, _family_key


def test_package_imports():
    import grains

    assert grains.__version__


def test_grains_csv_shape():
    df = pd.read_csv(config.GRAINS_CSV)
    expected = {
        "common_name",
        "species",
        "ncbi_taxon_id",
        "storage_protein_family",
        "celiac_status",
    }
    assert expected <= set(df.columns)
    assert df["celiac_status"].isin({"toxic", "safe", "intermediate"}).all()
    assert df["ncbi_taxon_id"].astype(str).str.isdigit().all()
    # the three gluten grains must be present and toxic
    toxic = set(df.loc[df["celiac_status"] == "toxic", "common_name"].str.lower())
    assert any("wheat" in t for t in toxic)
    assert "barley" in toxic and "rye" in toxic


def test_every_family_has_search_terms():
    df = pd.read_csv(config.GRAINS_CSV)
    for fam in df["storage_protein_family"]:
        key = _family_key(str(fam))
        assert key in config.STORAGE_PROTEIN_TERMS


def test_peptides_csv_shape():
    df = pd.read_csv(config.PEPTIDES_CSV)
    expected = {
        "epitope_id",
        "canonical_name",
        "restriction",
        "source_grain",
        "extended_sequence",
        "form",
        "category",
        "verified",
    }
    assert expected <= set(df.columns)
    assert df["restriction"].isin({"DQ2.5", "DQ8"}).all()
    assert df["form"].isin({"native", "deamidated"}).all()
    # extended sequences must be plausible peptide strings
    aa = set("ACDEFGHIKLMNPQRSTVWY")
    for seq in df["extended_sequence"].dropna():
        assert set(seq) <= aa, f"bad residues in {seq!r}"
        assert 8 <= len(seq) <= 40
    # at least the four immunodominant epitopes are marked verified
    dom = df[df["immunodominant"] == "yes"]
    assert (dom["verified"] == "yes").sum() >= 4


def test_core_is_a_substring_of_its_extended_sequence():
    """A row's own core must appear in the peptide that row says will be modeled.

    Regression for E09, which carried the deamidated core EGSFQPSQE alongside the *native*
    extended sequence QGSFQPSQQNPQ while labeled form=deamidated. Since extended_sequence is
    what gets modeled, that row would have modeled the native peptide under a deamidated
    label -- silently inverting the native-vs-deamidated comparison it exists to support.
    """
    df = pd.read_csv(config.PEPTIDES_CSV)
    cored = df[df["core_9mer"].notna() & (df["core_9mer"].astype(str) != "")]
    assert len(cored) >= 15, "expected most rows to declare a core"
    for _, row in cored.iterrows():
        assert row["core_9mer"] in row["extended_sequence"], (
            f"{row['epitope_id']} ({row['canonical_name']}): core "
            f"{row['core_9mer']!r} is not inside extended_sequence "
            f"{row['extended_sequence']!r}"
        )


def test_deamidated_rows_actually_carry_a_glutamate():
    """`form == deamidated` means at least one Q became E, so there must be an E.

    The converse of the E09 bug: a row claiming deamidation whose sequence has no glutamate
    isn't deamidated, whatever the column says.
    """
    df = pd.read_csv(config.PEPTIDES_CSV)
    for _, row in df[df["form"] == "deamidated"].iterrows():
        assert "E" in row["extended_sequence"], (
            f"{row['epitope_id']} is marked deamidated but its extended_sequence "
            f"{row['extended_sequence']!r} contains no glutamate"
        )


def test_native_deamidated_pairs_differ_only_by_q_to_e():
    """Each native_pair row must reduce to its partner by E->Q only, same length.

    This is what makes the pair a controlled comparison: if the two sequences differed by
    anything other than deamidation, a score difference between them wouldn't isolate the
    TG2 effect. Pairings come from the `notes` text ("Pairs with E0n").
    """
    df = pd.read_csv(config.PEPTIDES_CSV).set_index("epitope_id")
    pairs = {"N01": "E01", "N02": "E02", "N03": "E04", "C02": "C01"}
    for native_id, deamidated_id in pairs.items():
        native = df.loc[native_id, "extended_sequence"]
        deamidated = df.loc[deamidated_id, "extended_sequence"]
        assert len(native) == len(deamidated), f"{native_id}/{deamidated_id} length mismatch"
        differences = [
            (i + 1, n, d) for i, (n, d) in enumerate(zip(native, deamidated)) if n != d
        ]
        assert differences, f"{native_id} and {deamidated_id} are identical"
        for position, native_aa, deamidated_aa in differences:
            assert (native_aa, deamidated_aa) == ("Q", "E"), (
                f"{native_id}/{deamidated_id} differ at position {position} by "
                f"{native_aa}->{deamidated_aa}, which is not a TG2 deamidation"
            )


def test_build_query_syntax():
    q = build_query(4565, ["gliadin", "glutenin"])
    assert "taxonomy_id:4565" in q
    assert 'protein_name:"gliadin"' in q
    assert " OR " in q


@pytest.mark.parametrize(
    "family_field,expected",
    [
        ("gliadin/glutenin", "gliadin/glutenin"),
        ("hordein", "hordein"),
        ("prolamin", "prolamin"),
        ("something unknown", "prolamin"),
    ],
)
def test_family_key_mapping(family_field, expected):
    assert _family_key(family_field) == expected
