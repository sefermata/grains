# Checkpoint plan

Application due end of November 2026. Kickoff 2026-09-08. Checkpoints every two weeks;
each one ends with something you can show and talk about.

Time estimate assumes ~8–12 h/week; adjust as noted.

---

## Checkpoint 1 — Reference data & environment  (due 2026-09-22)

**Goal:** every input is downloaded, curated, cited, and sanity-checked.

- [x] `git init`, first commit. (2026-09-14 — used GitHub Desktop's bundled git instead of
      Xcode Command Line Tools; see README.md for why and how.)
- [x] `pip install -r requirements.txt`; `pytest` passes. (biopython 1.88 installed, 9/9 tests
      green as of 2026-09-16.)
- [x] `python -m grains.proteomes` — proteomes for all 16 grains in `data/raw/`,
      `data/processed/proteome_index.csv` written. Coverage gap for rye/teff/pearl-millet/
      amaranth confirmed real (public UniProt under-annotation) — noted, not hidden.
- [~] Verify **every** row of `data/reference/peptides.csv` against Sollid 2012/2020 + IEDB
      (see `data/reference/peptides_NOTES.md`). Add `iedb_id` and UniProt accession + range.
      **2026-09-16: name/core/restriction confirmed for all 17 rows against Sollid 2012**
      (caught and fixed two wrong PMIDs along the way). Still open: `iedb_id` column, UniProt
      accession+range column, the Sollid 2020 update, and a direct Tye-Din 2010 table check for
      N01–N03 (currently verified by derivation only).
- [x] Confirm the PDB target list (`data/reference/pdb_structures.md`); write
      `data/processed/targets.csv`. Download the chosen validation structure (`4OZF` or `1S9V`).
      **2026-09-21: confirmed all 8 candidate PDB IDs against rcsb.org.** Two corrections to
      the from-memory list: `4OZF/4OZG/4OZH/4OZI` are all ternary (TCR-bound), not binary as
      guessed, so `1S9V` — not `4OZF` — is the Checkpoint 2 validation structure; and `4D0N`
      isn't a gluten/DQ structure at all (wrong ID) while `2NNA` turned out to genuinely be a
      DQ8+gluten binary structure (previously guessed as an unrelated insulin reference).
      Bonus: by exact sequence match, `1S9V`'s peptide = panel epitope E01 and `2NNA`'s = E09,
      so validation directly checks panel peptides, not off-panel ones. `1S9V.pdb` and
      `2NNA.pdb` downloaded to `data/raw/pdb/`.
- [x] Render one crystal structure (ChimeraX / PyMOL / nglview) showing a gluten peptide in
      the DQ2.5 groove — save to `results/figures/`. **2026-09-21: no ChimeraX/PyMOL install**
      (disk-constrained, see README) — used a Biopython CA-backbone-trace rendering instead
      (`src/grains/render_structure.py` → `results/figures/1S9V_dq25_groove.png`). Good enough
      to show the peptide nestled in the groove; a nicer cartoon render is a Checkpoint 5 polish
      item if disk space frees up.

**Present:** README + `grains.csv` + verified `peptides.csv` + the rendered structure, with a
2–3 sentence statement of the question and the validation logic.

---

## Checkpoint 2 — Pipeline validation  (due 2026-10-06)

**Goal:** show the modeling approach can reproduce a known complex.

- [x] `fetch_hla.py` — extract DQA1*05:01 + DQB1*02:01 groove+Ig domains → `hla_dq25.fasta`.
      **2026-09-21:** reads chains straight from `1S9V.pdb`'s `SEQRES` records (the crystallized
      construct already *is* the groove+Ig ectodomain, no domain-boundary guessing needed) —
      not from full-length UniProt entries. Regression-tested against a real 7-residue disordered
      loop gap in chain B that a naive CA-atom walk would silently splice across.
- [ ] Submit the validation peptide + HLA chains to AlphaFold Server; retrieve the model.
- [ ] Superpose model on the crystal structure (Biopython `Superimposer` or ChimeraX
      `matchmaker`); report backbone RMSD for (a) whole complex, (b) peptide core.
- [ ] Try the CPU-only homology route as a cross-check (PANDORA) on the same peptide.
- [ ] Write `src/grains/modeling.py` wrapping "sequence(s) in → model file out" for both routes.

**AlphaFold Server Terms-of-Service constraint (found 2026-09-21, see Checkpoint 3):**
AlphaFold Server output may be compared to the crystal structure by RMSD here — that's a pure
geometric check, not a binding prediction, so it's fine. Its output must **never** be piped into
PRODIGY or any other binding/affinity scorer (that's a prohibited use — see Checkpoint 3). Keep
this checkpoint's AlphaFold Server usage scoped to exactly that: one validation structure, one
RMSD comparison, nothing downstream of it feeds a binding score.

**Present:** overlay figure (model vs crystal) + RMSD numbers + a short method paragraph.
*If AlphaFold Server access is slow:* present the PANDORA result and the submission plan.

---

## Checkpoint 3 — Full panel + scoring  (due 2026-10-20)

**Goal:** every panel peptide modeled and scored.

**Constraint, locked in 2026-09-21 (AlphaFold Server Prohibited Use Policy, clause 2):**
AlphaFold Server output must not be used "in connection with any automated system that predicts
the binding or interaction of the protein with ligands or peptides" — PRODIGY is exactly that
kind of system (protein–peptide binding-affinity prediction), same category as the policy's own
named examples (Glide, AutoDock). So: **every one of the ~40-50 panel complexes that will be
PRODIGY-scored must be modeled with PANDORA (CPU-only, MODELLER-based — a separate, unrestricted
pipeline), never with AlphaFold Server.** AlphaFold Server stays confined to the one-time
Checkpoint 2 validation structure, whose only downstream use is an RMSD comparison — never a
PRODIGY input. If PANDORA's accuracy on the validation peptide (Checkpoint 2) is too poor to
trust for the full panel, that's a real limitation to state plainly in the report — not a reason
to reach for AlphaFold Server here instead.

- [ ] `build_controls.py` — length- and P/Q-matched control peptides from zein / kafirin /
      orzenin / avenin + a non-epitope gluten stretch → `data/processed/controls.csv`.
- [ ] Model toxic panel + controls + native pairs (~40–50 complexes) **with PANDORA**.
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
- [ ] Replace `results/figures/1S9V_dq25_groove.png` (Checkpoint 1's quick Biopython
      CA-backbone-trace render) with a real ChimeraX or PyMOL cartoon: helices as ribbons,
      beta-sheet as arrows, peptide as sticks colored by atom, viewed straight down into the
      groove. The current figure only shows the peptide is roughly in the right place; it
      can't show anchor-pocket chemistry (the deamidated Glu side chains dropping into the
      DQ2.5 pockets), which is the actual story worth a headline figure. Needs disk space
      freed up first (no Xcode CLT / low disk currently blocks a ChimeraX/PyMOL install).
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
