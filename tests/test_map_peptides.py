"""Tests for locating panel peptides in real UniProt sequences.

The behavior worth pinning down is the deamidation-tolerant match: a deamidated peptide is
not a literal substring of any native gene, so a literal search would report every deamidated
epitope as absent. Allowing peptide E to match gene Q -- and only in that direction -- both
finds them and identifies the deamidation sites.
"""

from __future__ import annotations

import pytest

from grains import map_peptides
from grains.map_peptides import Protein

pytestmark = pytest.mark.skipif(
    not (map_peptides.config.RAW / "bread-wheat__gliadin-glutenin.fasta").exists(),
    reason="proteome fastas not downloaded (data/raw is git-ignored; run grains.proteomes)",
)


def _protein(sequence: str, grain: str = "wheat") -> Protein:
    return Protein("P00000", "TEST_PROT", "test protein", sequence, grain)


def test_literal_match_reports_no_deamidation():
    protein = _protein("AAAPFPQPELPYAAA")
    start, flexed = map_peptides.find_in_protein("PFPQPELPY", protein)
    assert start == 4
    assert flexed == []


def test_deamidated_peptide_matches_native_gene_and_names_the_site():
    """E01's core against the native gene: the E at core position 6 comes from a gene Q."""
    protein = _protein("AAAPFPQPQLPYAAA")  # native: Q where the epitope has E
    start, flexed = map_peptides.find_in_protein("PFPQPELPY", protein)
    assert start == 4
    assert flexed == [6]


def test_q_is_not_allowed_to_match_gene_e():
    """Deamidation runs Q->E in the peptide relative to the gene, never the reverse.

    Permitting the reverse would let a native peptide 'match' a gene that encodes glutamate,
    inventing a deamidation site that chemistry doesn't support.
    """
    protein = _protein("AAAPFPQPELPYAAA")  # gene has E
    assert map_peptides.find_in_protein("PFPQPQLPY", protein) is None


def test_non_deamidation_mismatch_is_rejected():
    protein = _protein("AAAPFPQPELPWAAA")  # final Y -> W
    assert map_peptides.find_in_protein("PFPQPELPY", protein) is None


def test_pattern_only_flexes_glutamate():
    pattern = map_peptides.deamidation_tolerant_pattern("PEQ")
    assert pattern.search("PEQ")
    assert pattern.search("PQQ")      # E may be Q
    assert not pattern.search("PEE")  # Q may not be E


def test_known_epitopes_land_where_expected():
    """E01 and its native partner N01 must map to the same gene window.

    A native/deamidated pair that mapped to different places wouldn't be a controlled pair,
    and the TG2 comparison built on it would be comparing two different loci.
    """
    rows = {row["epitope_id"]: row for row in map_peptides.run()}

    e01, n01 = rows["E01"], rows["N01"]
    assert e01["uniprot_accession"] == n01["uniprot_accession"]
    assert (e01["start"], e01["end"]) == (n01["start"], n01["end"])
    assert e01["match_type"] == "deamidated"
    assert n01["match_type"] == "exact"
    assert e01["n_deamidation_sites"] == 1


def test_33mer_has_three_deamidation_sites():
    """C01 carries three embedded glia-a1a-type cores, so three Q->E relative to C02."""
    rows = {row["epitope_id"]: row for row in map_peptides.run()}
    assert rows["C01"]["n_deamidation_sites"] == 3
    assert rows["C02"]["match_type"] == "exact"
    assert (rows["C01"]["start"], rows["C01"]["end"]) == (
        rows["C02"]["start"], rows["C02"]["end"]
    )


def test_every_deamidated_row_needs_at_least_one_flex():
    """A row marked deamidated that matches a gene literally is mislabeled or has the
    wrong sequence -- this is what caught E09 carrying a native extended_sequence."""
    for row in map_peptides.run():
        if row["match_type"] == "exact":
            assert row["form"] == "native", (
                f"{row['epitope_id']} matched a gene literally (no Q->E needed) but is "
                f"marked form={row['form']!r} -- a deamidated peptide must differ from the "
                "gene by at least one Q->E"
            )
