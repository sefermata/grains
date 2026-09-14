# grains — structural screening of cereal peptides against HLA-DQ2.5

**Question:** Do grain-derived peptides bound to HLA-DQ2.5 explain which cereals are toxic
in celiac disease? We model known immunogenic gluten peptides and length-/composition-matched
control peptides from *safe* grains inside the HLA-DQ2.5 peptide-binding groove, score the
predicted binding, and test whether structure-based scores separate toxic from safe grains —
benchmarked against known biology and a sequence-based predictor.

## Background (plain language)

In celiac disease the immune system reacts to **gluten**, a family of seed-storage proteins in
wheat, barley, and rye. The disease-driving events are molecular and well characterized:

1. Gluten proteins are proline/glutamine-rich and resist digestion, so short peptides
   (~9–20 residues) survive intact in the gut.
2. The enzyme **transglutaminase 2 (TG2)** converts specific glutamine (Q) residues to
   glutamate (E) — "deamidation".
3. Deamidated peptides bind tightly in the peptide-binding groove of **HLA-DQ2.5**
   (or HLA-DQ8), a class II MHC molecule that displays peptides to CD4+ T cells.
4. T cells recognize the displayed peptide as foreign and drive the inflammatory response
   that damages the small intestine.

Grains that are safe in celiac (rice, maize, sorghum, millets, and — with caveats — oats)
either lack these peptides or carry very few. Because the toxic/safe status of every common
grain is already known, this project has a **built-in ground truth** to validate against.

The molecular recognition step we model is step 3: **peptide ↔ HLA-DQ2.5 groove**. There is no
"gluten receptor"; HLA-DQ2.5 is the relevant target, and the PDB contains crystal structures of
it bound to gluten peptides that serve as gold-standard references.

## Approach

| Phase | Work | Presentable as |
|------|------|----------------|
| 1 | Reference data: grain list + taxonomy, curated toxic-peptide panel, safe-grain proteomes from UniProt, target HLA-DQ2.5 structures from the PDB | README + cited data tables + a rendered DQ2.5–gluten crystal structure |
| 2 | Pipeline validation: rebuild a peptide from a known crystal structure, compare to the real complex (RMSD) | "reproduce a known DQ2.5 complex to *X* Å" + overlay figure |
| 3 | Model the full peptide panel (AlphaFold Server), compute binding scores (interface contacts, buried surface area, PRODIGY affinity, model confidence) | score table + toxic-vs-safe comparison plot |
| 4 | Sequence baseline (IEDB MHC-II API), statistics on toxic vs safe, effect of TG2 deamidation | figures + stats + interpretation |
| 5 | Writeup (2–3 pp), clean repo, end-to-end notebook, limitations | finished portfolio piece |

See [`docs/checkpoints.md`](docs/checkpoints.md) for the dated checkpoint plan and task lists.

## Repository layout

```
data/
  reference/     curated, version-controlled inputs (grains.csv, peptides.csv, pdb_structures.md)
  raw/           downloaded proteomes / structures (not version-controlled)
  processed/     derived tables
src/grains/      pipeline code
  config.py      paths, UniProt settings, search keywords, HLA accessions
  proteomes.py   fetch seed-storage proteins from UniProt
notebooks/       exploratory analysis
results/         figures and score tables
tests/           smoke tests for the reference data and code
docs/            checkpoint plan, notes
```

## Setup

```bash
# Anaconda Python 3.13 already provides pandas/numpy/matplotlib/scipy/sklearn/seaborn.
pip install -r requirements.txt        # adds biopython (+ scoring tools later)

# fetch safe/toxic grain seed-storage proteomes from UniProt
python -m grains.proteomes
```

## Data sources

- **Gluten T-cell epitopes:** Sollid LM *et al.*, "Nomenclature and listing of celiac disease
  relevant gluten epitopes recognized by CD4+ T cells", *Immunogenetics* 2012 (PMID 22371225)
  and the 2020 update (PMID 31735986). Cross-checked against the IEDB (iedb.org).
- **Proteomes:** UniProtKB REST API (rest.uniprot.org).
- **Structures:** RCSB PDB (rcsb.org) — see `data/reference/pdb_structures.md`.
- **Sequence-based binding baseline:** IEDB MHC class II binding prediction API.
- **Structure prediction:** AlphaFold Server (alphafoldserver.com).
- **Affinity scoring:** PRODIGY (bianca-cimr / haddocking).

## Limitations & disclaimer

This is a **research and educational** project. It models one step (peptide–HLA-DQ2.5 binding)
of a multi-step disease mechanism and ignores TG2 substrate specificity in detail, T-cell
receptor recognition, innate immune contributions (e.g. the p31–43 peptide), amylase-trypsin
inhibitors, dose, and real-world cross-contamination. In-silico binding scores are not
measurements and are not a basis for any dietary or medical decision. Nothing here is medical
advice.
