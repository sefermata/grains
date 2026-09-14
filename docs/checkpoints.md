# Checkpoint plan

Application due end of November 2026. Kickoff 2026-09-08. Checkpoints every two weeks;
each one ends with something you can show and talk about.

Time estimate assumes ~8–12 h/week; adjust as noted.

---

## Checkpoint 1 — Reference data & environment  (due 2026-09-22)

**Goal:** every input is downloaded, curated, cited, and sanity-checked.

- [ ] Install Xcode Command Line Tools (`xcode-select --install`), then `git init`, first commit.
- [ ] `pip install -r requirements.txt`; `pytest` passes.
- [ ] `python -m grains.proteomes` — proteomes for all 16 grains in `data/raw/`,
      `data/processed/proteome_index.csv` written. Note which grains returned few/no hits.
- [ ] Verify **every** row of `data/reference/peptides.csv` against Sollid 2012/2020 + IEDB
      (see `data/reference/peptides_NOTES.md`). Add `iedb_id` and UniProt accession + range.
- [ ] Confirm the PDB target list (`data/reference/pdb_structures.md`); write
      `data/processed/targets.csv`. Download the chosen validation structure (`4OZF` or `1S9V`).
- [ ] Render one crystal structure (ChimeraX / PyMOL / nglview) showing a gluten peptide in
      the DQ2.5 groove — save to `results/figures/`.

**Present:** README + `grains.csv` + verified `peptides.csv` + the rendered structure, with a
2–3 sentence statement of the question and the validation logic.

---

## Checkpoint 2 — Pipeline validation  (due 2026-10-06)

**Goal:** show the modeling approach can reproduce a known complex.

- [ ] `fetch_hla.py` — extract DQA1*05:01 + DQB1*02:01 groove+Ig domains → `hla_dq25.fasta`.
- [ ] Submit the validation peptide + HLA chains to AlphaFold Server; retrieve the model.
- [ ] Superpose model on the crystal structure (Biopython `Superimposer` or ChimeraX
      `matchmaker`); report backbone RMSD for (a) whole complex, (b) peptide core.
- [ ] Try the CPU-only homology route as a cross-check (PANDORA) on the same peptide.
- [ ] Write `src/grains/modeling.py` wrapping "sequence(s) in → model file out" for both routes.

**Present:** overlay figure (model vs crystal) + RMSD numbers + a short method paragraph.
*If AlphaFold Server access is slow:* present the PANDORA result and the submission plan.

---

## Checkpoint 3 — Full panel + scoring  (due 2026-10-20)

**Goal:** every panel peptide modeled and scored.

- [ ] `build_controls.py` — length- and P/Q-matched control peptides from zein / kafirin /
      orzenin / avenin + a non-epitope gluten stretch → `data/processed/controls.csv`.
- [ ] Model toxic panel + controls + native pairs (~40–50 complexes).
- [ ] `src/grains/scoring.py` — per-complex: interface contact count, buried surface area
      (freesasa), PRODIGY ΔG / Kd, interface PAE + peptide pLDDT, anchor-pocket occupancy.
- [ ] `data/processed/scores.csv`; first comparison plot (toxic vs safe vs intermediate).

**Present:** score table + a boxplot/strip plot of one or two scores by celiac status.

---

## Checkpoint 4 — Analysis & baseline  (due 2026-11-03)

**Goal:** turn scores into a defensible answer.

- [ ] Sequence baseline: IEDB MHC-II API (`netmhciipan` method, alleles HLA-DQA1*05:01/
      DQB1*02:01 and DQ8) for every panel peptide → `data/processed/seq_baseline.csv`.
- [ ] `src/grains/deamidation.py` — TG2 QxP-type rules; check the `form` pairs are consistent.
- [ ] Stats: Mann–Whitney U (toxic vs safe) per score; effect sizes; deamidated vs native
      (paired). ROC / AUC using celiac status as the label.
- [ ] Compare structure scores vs the sequence baseline — where does structure add signal?

**Present:** 2–3 figures with statistics + a paragraph interpreting them.

---

## Checkpoint 5 — Writeup & polish  (due 2026-11-17)

**Goal:** a portfolio-ready artifact.

- [ ] `notebooks/analysis.ipynb` runs top to bottom from `data/processed/`.
- [ ] `report.md` (2–3 pp): question, methods, results, limitations, "what I'd do next".
- [ ] README updated with results summary + a headline figure.
- [ ] Repo tidy: docstrings, `pytest` green, `requirements.txt` complete, LICENSE.
- [ ] Optional stretch: a second HLA allele (DQ8), or ancient-wheat (einkorn/emmer) prolamins.

**Present:** the repo + report as the portfolio piece.

**2026-11-17 → 11-28:** buffer for slippage / stretch goals / application text.

---

## Risk log

| Risk | Mitigation |
|------|-----------|
| AlphaFold Server rate limits / queue | PANDORA (CPU) is the fallback modeling route; start submissions early |
| Sparse UniProt coverage for teff / millets | Drop to a smaller grain set; note the data gap explicitly in the report |
| Peptide sequence errors propagate | Checkpoint 1 verification gate; `verified == yes` required before modeling |
| Scope creep | TCR step, MD refinement, and extra alleles are explicitly out of scope for v1 |
