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
