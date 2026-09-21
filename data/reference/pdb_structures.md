# Target structures — HLA-DQ2.5 (and HLA-DQ8) with gluten peptides

HLA-DQ2.5 = the heterodimer of **DQA1\*05:01** (alpha chain) + **DQB1\*02:01** (beta chain).
HLA-DQ8 = DQA1\*03:01 + DQB1\*03:02.

**Confirmed against rcsb.org 2026-09-21** (Checkpoint 1). Full list with resolution,
has_TCR, and exact peptide sequence is in `data/processed/targets.csv`. Corrections vs.
the from-memory list below:
- `4OZF`/`4OZG`/`4OZH`/`4OZI` are **all ternary** (TCR-bound) — none is the "binary,
  high-res, single epitope" structure this doc originally guessed for 4OZF. Checked by
  chain composition (alpha, beta, TCR-alpha, TCR-beta, peptide = 5 entities each).
- `4OZI`'s peptide sequence (`QPFPQPELPY`) matches the **glia-α1a (E01)** core, not
  glia-α2 as RCSB's own structure-title text implies — verified by direct sequence
  comparison against `data/reference/peptides.csv`, not by trusting the summary text.
- `4D0N` is **not a gluten/DQ structure at all** — it's an unrelated AKAP13/RhoA GEF
  complex. Wrong ID, dropped from the list entirely.
- `2NNA` **is** a real DQ8 + gluten peptide binary structure (this doc previously
  guessed it was an unrelated insulin-B reference — wrong). Its peptide
  (`QQYPSGEGSFQPSQENPQ`) contains the **E09 (DQ8-glia-a1) core exactly**
  (`EGSFQPSQE`), so it doubles as a direct match to a panel epitope.

**Checkpoint 2 validation structure: `1S9V`.** Binary (no TCR), 2.2 Å, and its peptide
(`LQPFPQPELPY`) matches the **E01 (DQ2.5-glia-a1a)** core exactly — it directly validates
a panel peptide rather than an off-panel one. `4OZF` is disqualified as "first choice"
since it's ternary; there is no binary, non-TCR DQ2.5+single-epitope structure among the
Petersen 2014 series (4OZF–4OZI) — they're all TCR-bound. `2NNA` is the equivalent binary
validation structure if/when a DQ8 cross-check is wanted (matches E09).

Rendered: `results/figures/1S9V_dq25_groove.png` (CA backbone trace, peptide in red,
nestled between the alpha/beta chains — see `src/grains/render_structure.py`).

## HLA-DQ2.5 + gluten peptide (primary references)

| PDB | Peptide (deposited) | Epitope match | has_TCR | Resolution | Reference |
|-----|---------------------|----------------|---------|------------|-----------|
| `1S9V` | `LQPFPQPELPY` | E01 (glia-α1a) | no | 2.2 Å | Kim et al. 2004, PNAS (PMID 15020763) |
| `4OZF` | `APQPELPYPQP` (+PGS tag) | E02 (glia-α2) | yes (JR5.1) | 2.7 Å | Petersen et al. 2014, NSMB (PMID 24777060) |
| `4OZG` | `APQPELPYPQP` (+PGS tag) | E02 (glia-α2) | yes (d2) | 3.0 Å | Petersen et al. 2014 |
| `4OZH` | `APQPELPYPQP` (+PGS tag) | E02 (glia-α2) | yes (s16) | 2.8 Å | Petersen et al. 2014 |
| `4OZI` | `QPFPQPELPY` (+PGS tag) | E01 (glia-α1a) | yes (s2) | 3.2 Å | Petersen et al. 2014 |

## HLA-DQ8 + gluten peptide

| PDB | Peptide (deposited) | Epitope match | has_TCR | Resolution | Reference |
|-----|---------------------|----------------|---------|------------|-----------|
| `2NNA` | `QQYPSGEGSFQPSQENPQ` | E09 (DQ8-glia-a1) | no | 2.1 Å | Henderson et al. 2007, Immunity (PMID 17629515) |
| `4GG6` | `QQYPSGEGSFQPSQENPQ` | E09 (DQ8-glia-a1) | yes (SP3.4) | 3.2 Å | Broughton et al. 2012, Immunity (PMID 23063329) |

## How these are used

- **Checkpoint 2 (validation):** use `1S9V` — the only binary (peptide + MHC only, no TCR)
  DQ2.5 structure on this list (see correction above: 4OZF/G/H/I all turned out to be
  ternary). Extract the peptide sequence and the two MHC chains, rebuild the complex with
  the modeling pipeline, superpose on the crystal structure, report backbone RMSD over
  (a) the whole complex and (b) the peptide core only.
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
