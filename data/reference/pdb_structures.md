# Target structures — HLA-DQ2.5 (and HLA-DQ8) with gluten peptides

HLA-DQ2.5 = the heterodimer of **DQA1\*05:01** (alpha chain) + **DQB1\*02:01** (beta chain).
HLA-DQ8 = DQA1\*03:01 + DQB1\*03:02.

All PDB IDs below are from memory and **must be confirmed on rcsb.org during Checkpoint 1**
(check: does the entry contain a gluten-derived peptide? is it DQ2.5 or DQ8? resolution?).
Record the confirmed list in `data/processed/targets.csv` with columns
`pdb_id, allele, peptide_seq, peptide_epitope_id, resolution_A, has_TCR, reference_pmid`.

## HLA-DQ2.5 + gluten peptide (primary references)

| PDB (verify) | Peptide | Notes | Reference |
|-------------|---------|-------|-----------|
| `1S9V` | α-II gliadin (deamidated) | First DQ2.5–gliadin structure | Kim et al. 2004, PNAS (PMID 15020763) |
| `4OZF` | DQ2.5-glia-α1a (deamidated) | High-res DQ2.5 + single epitope | Petersen et al. 2014, Nat Struct Mol Biol (PMID 24777060) |
| `4OZG` | DQ2.5-glia-α2 (deamidated) | Companion to 4OZF | Petersen et al. 2014 |
| `4OZH` / `4OZI` | α1a or α2 + TCR | Ternary complex (adds a TCR chain — heavier to model) | Petersen et al. 2014 |
| `4D0N` | deamidated gliadin + TCR | Ternary complex | Petersen et al. 2016 (verify) |

## HLA-DQ8 + gluten peptide

| PDB (verify) | Peptide | Notes | Reference |
|-------------|---------|-------|-----------|
| `2NNA` | insulin B (NOT gluten) | apo-ish DQ8 reference only | — |
| `4GG6` | deamidated gliadin + TCR | DQ8 ternary complex (verify id/allele) | Broughton et al. 2012 |

## How these are used

- **Checkpoint 2 (validation):** pick one *binary* (peptide + MHC only, no TCR) structure —
  `4OZF` is the first choice, `1S9V` the backup. Extract the peptide sequence and the two MHC
  chains, rebuild the complex with the modeling pipeline, superpose on the crystal structure,
  report backbone RMSD over (a) the whole complex and (b) the peptide core only.
- **Checkpoint 3 (panel):** use the DQ2.5 alpha+beta chains from the chosen reference as the
  fixed receptor; thread every panel peptide through the groove.
- TCR-containing structures are **not** modeled — this project stops at peptide–MHC.

## Chain sequences

Do not hard-code these. Fetch at Checkpoint 2:
- DQA1\*05:01 mature ectodomain — from the chosen PDB entry, or UniProt (HLA-DQ alpha chain).
- DQB1\*02:01 mature ectodomain — from the chosen PDB entry, or UniProt `P01920` region.
- Keep only the alpha-1/alpha-2 (or beta-1/beta-2) domains that form the groove + Ig fold;
  drop signal peptide, transmembrane, and cytoplasmic segments. `fetch_hla.py` (Checkpoint 2)
  will do this and write `data/processed/hla_dq25.fasta`.
