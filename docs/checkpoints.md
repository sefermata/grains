# Checkpoint plan

Application due end of November 2026. Kickoff 2026-09-08. Each checkpoint ends with something
you can show and talk about.

New in the project? Read [`background.md`](background.md) first — the biology, the question, and
what each pipeline stage is for.

Time estimate assumes ~8–12 h/week; adjust as noted.

**Timeline pulled in 2026-10-06.** Original plan ran Checkpoint 5 to 2026-11-17 with eleven
days of buffer. Revised: **everything substantive is done by 2026-11-07**, and 11-07 → 11-28 is
finishing touches and fixes. Checkpoints 1 and 2 are unchanged (done); Checkpoint 3 keeps its
10-20 date because it is the heaviest compute lift and already scheduled tight. The compression
lands on Checkpoints 4 and 5:

| Checkpoint | Was | Now |
|---|---|---|
| 3 — full panel + scoring | 2026-10-20 | **2026-10-20** (unchanged) |
| 4 — analysis & baseline | 2026-11-03 | **2026-10-30** |
| 5 — writeup & polish | 2026-11-17 | **2026-11-07** |
| finishing touches / fixes | 11-17 → 11-28 (11 days) | **11-07 → 11-28 (3 weeks)** |

*Why this is an improvement, not just a squeeze:* a report written once against a deadline is
worse than a draft written early and revised for three weeks. Moving Checkpoint 5 to 11-07
converts buffer-for-slippage into revision time, which is where a portfolio piece actually
gets better. The cost is that **Checkpoint 3 now has no slack** — it was the one unknown-duration
stage (PANDORA on ~40–50 complexes, disk-constrained), so it gets an explicit go/no-go gate
below rather than being allowed to drift.

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
      (caught and fixed two wrong PMIDs along the way).

      **2026-10-05: steps 2, 3 and 4 done; the checklist is closed bar one naming question.**
      `iedb_id` (step 2, `grains.iedb`, 12/15 cored rows) and `uniprot_accession` +
      `uniprot_range` (step 3, `grains.map_peptides`) are now columns in `peptides.csv`, and
      the mapping recovered the deamidation positions too (step 4). All 17 rows are
      gene-faithful: 13 deamidated rows each need ≥1 Q→E against the gene, 4 native rows match
      literally — that split is itself the check.

      Four data errors came out of it, all of the silent kind, and all fixed:
      **three more wrong PMIDs** (including `22371225`, cited by every epitope row, which
      actually points at an organic-chemistry paper → `22322673`); **E09 carrying a native
      `extended_sequence`** under a `deamidated` label, which would have inverted the TG2
      comparison since that field is the modeling input; and **five rows with fabricated
      flanks** (E05, E06, E08, E11, E12) whose verified cores sit in real genes but whose
      surrounding residues match no gene at all. Full account in `peptides_NOTES.md`.

      **Still open (one item):** the E10 `DQ8-glia-γ1a` vs `g1` name-suffix question, which
      needs the Sollid 2020 text.
      **Closed as limitations, not tasks:** Tye-Din 2010 and Sollid 2020 are confirmed
      paywalled (no PMC record), so N01–N03 stay derivation-verified — though the UniProt
      mapping independently showed each sits at the *same gene locus* as its deamidated
      partner differing only by Q→E, which is stronger evidence than the original derivation.
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
- [x] Submit the validation peptide + HLA chains to AlphaFold Server; retrieve the model.
      **2026-10-05: done — the model reproduces the crystal structure.** 5 models returned,
      all pTM/ipTM 0.91. Core CA RMSD to 1S9V: **0.59–0.73 Å**, i.e. *better* than two
      independent crystal structures of the same complex agree with each other (0.91–1.16 Å),
      sitting between that ceiling and the 0.24 Å within-crystal refinement floor. Groove fit
      0.69–0.89 Å. **The pipeline reproduces a known complex — Checkpoint 2's question is
      answered yes.**

      **Not circular, and checked rather than assumed:** 1S9V is *not* among the templates
      AlphaFold used (6px6, 4ozi, 4d8p, 6dig, 6u3m, 3pl6, 1uvq, 1iak, 4grl). **4OZI is**, for
      the peptide chain — a caveat worth stating, so it was tested rather than waved at. Had
      the model simply recalled its template it would sit at 4OZI's own distance from 1S9V
      (~1.1 Å core) and be much closer to 4OZI than to 1S9V. Measured: it is roughly
      *equidistant* (core 0.59–0.73 Å to 1S9V, 0.48–0.73 Å to 4OZI) while the two crystals are
      0.91–1.16 Å apart from each other — a consensus prediction, not template recall. It also
      fits the binary 1S9V groove (0.69–0.89 Å) far better than the TCR-bound 4OZI groove
      (1.53–1.67 Å), i.e. it reproduced the binary conformation rather than the ternary one it
      was templated from.

      Reproduce with `PYTHONPATH=src python -c "from grains.superpose import validate_models;
      validate_models()"` → `data/processed/rmsd_validation.csv`.
      Submission prep (same day): **submission file written** — `data/processed/alphafold_validation_job.json`
      (3 chains: DQA1\*05:01 193 aa + DQB1\*02:01 198 aa + peptide `LQPFPQPELPY`, the
      sequence as deposited in 1S9V chain C so the RMSD needs no sequence caveat).
      Format confirmed against AlphaFold Server's own JSON spec (`alphafoldserver` dialect,
      version 1). No PTM declared for the deamidated residue — TG2 turns Gln into Glu, and
      Glu is a standard amino acid, so it's just the `E` in the sequence.
      **Blocked on a manual upload**: AlphaFold Server has no public API, so the job has to
      go through the web UI by hand (bulk-upload/import JSON), then the mmCIF comes back to
      `data/raw/models/`.
- [x] Superpose model on the crystal structure (Biopython `Superimposer` or ChimeraX
      `matchmaker`); report backbone RMSD for (a) whole complex, (b) peptide core.
      **2026-10-05: `src/grains/superpose.py` written and the crystal-vs-crystal baseline
      measured** (`data/processed/rmsd_baseline.csv`). An RMSD is uninterpretable on its own
      — "is 1.5 Å good?" has no answer without knowing how well crystal structures of the
      same complex agree with *each other* — so the noise floor was measured first, from
      structures alone, before any model exists:

      | Comparison | Groove fit | Peptide | E01 core |
      |---|---|---|---|
      | 1S9V copy 1 vs 2 (same crystal, 2.2 Å) | 0.40 Å | 0.30 Å | **0.24 Å** |
      | 4OZI copy 1 vs 2 (same crystal, 3.2 Å) | 0.75 Å | 0.64 Å | **0.57 Å** |
      | 1S9V vs 4OZI (independent depositions) | 1.40 Å | 1.14 Å | **1.16 Å** |
      | 1S9V vs 4OZI (second copies) | 0.92 Å | 0.88 Å | **0.91 Å** |

      So: **a model within ~1.2 Å on the E01 core has matched experimental reproducibility**;
      ~0.25 Å is the refinement-noise floor and nothing can be expected to beat it. Beyond
      ~2 Å the model disagrees with the crystal more than independent crystals of the same
      complex disagree with each other.

      Protocol: fit on the MHC groove (α+β CA atoms), then measure the peptide *where it
      lands*. Never re-fit on the peptide — that reports the peptide's internal shape while
      discarding whether it sits in the same place in the groove, which is the whole question.
      Residue correspondence is by sequence alignment, never residue number (the α chains
      carry insertion codes; the peptide is numbered 1–11 in 1S9V but 2–12 in 4OZI).

      **Bonus finding, relevant to scope:** the peptide conformation barely changes between
      the binary complex (1S9V, no TCR) and the TCR-bound one (4OZI) — ~1 Å on the core,
      within the cross-deposition spread. That's direct evidence for the v1 scope decision to
      model peptide–MHC without the TCR: the TCR is not reshaping the bound peptide.
- [ ] Try the CPU-only homology route as a cross-check (PANDORA) on the same peptide.
      **Blocked 2026-10-05:** neither PANDORA nor MODELLER is installed. `csb-pandora` 0.9 is
      on PyPI and `modeller` 10.8 has a py313 build on the `salilab` conda channel, so the
      install path is clear, but **MODELLER needs a free academic license key** (register at
      <https://salilab.org/modeller/registration.html>) — a manual step. This is the critical
      path: Checkpoint 3 locks PANDORA in as the *only* route for complexes that get scored.
      `grains.modeling.pandora_available()` reports readiness and the exact install commands.
- [x] Write `src/grains/modeling.py` wrapping "sequence(s) in → model file out" for both routes.
      **2026-10-05:** AlphaFold-Server half is done and tested (builds the submission JSON).
      PANDORA half is an availability check + install instructions until the license key
      exists; the model-running call is deliberately not written blind against an API that
      can't be executed here yet.
      The module also turns the ToS constraint below into a **runtime guard**,
      `assert_scorable(route)`: Checkpoint 3's scoring step calls it per complex, and it
      refuses any non-PANDORA route (default-deny, so an unvetted new route isn't scorable
      just by being new). The failure mode it prevents is invisible — a PRODIGY number
      computed from an AlphaFold model looks exactly like a valid one.

**AlphaFold Server Terms-of-Service constraint (found 2026-09-21, see Checkpoint 3):**
AlphaFold Server output may be compared to the crystal structure by RMSD here — that's a pure
geometric check, not a binding prediction, so it's fine. Its output must **never** be piped into
PRODIGY or any other binding/affinity scorer (that's a prohibited use — see Checkpoint 3). Keep
this checkpoint's AlphaFold Server usage scoped to exactly that: one validation structure, one
RMSD comparison, nothing downstream of it feeds a binding score.

**Present:** overlay figure (model vs crystal) + RMSD numbers + a short method paragraph.
**Figure done 2026-10-05:** `results/figures/af_vs_1S9V_overlay.png` (+ `_dark.png`), from
`grains.render_overlay`. Two panels: the superposed peptide backbones, and per-residue CA
deviation against the two measured reference levels. The per-residue panel carries the real
result — deviation is *lowest in the 9-mer core* (0.13–0.40 Å, at or below the within-crystal
floor) and highest at the two flanking termini (0.83, 1.18 Å). That is exactly the MHC class II
expectation: the core is clamped in the groove while the flanks extend into solvent, so the
model is most accurate precisely where the binding happens. The reverse pattern would have been
a warning.
*If AlphaFold Server access is slow:* present the PANDORA result and the submission plan.

**Status 2026-10-05:** the measuring apparatus and its interpretive baseline are done and
tested; what's missing is the two models to put through it, each blocked on a manual step
(AlphaFold Server upload; MODELLER license key). The crystal-vs-crystal baseline is
presentable on its own — it's the number that makes any later RMSD mean something.

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

**Carry-forward from Checkpoint 1 (2026-10-05), both affecting the statistics:**
- **E03 and E06 are not independent observations.** DQ2.5-glia-ω1 (wheat ω-gliadin,
  `PFPQPEQPF`) and DQ2.5-sec-2 (rye secalin, `PQPEQPFPQ`) overlap by **7 of 9 residues** —
  two offset registers of the same repetitive prolamin motif, and no IEDB entry separates
  them. They will model and score almost identically, so counting both toward a toxic-vs-safe
  comparison overstates the sample. Treat them as one observation, or report with and without
  E06 and say which was used. (It is also part of *why* rye is toxic — worth a line in the
  writeup rather than only a caveat.)
- **Peptide length is a confound.** Panel `extended_sequence` lengths run 12–18 aa, and a
  longer peptide buries more surface whatever its chemistry. Match controls per peptide
  rather than to a panel average, and report length as a covariate or use per-residue
  normalized scores — never compare raw buried surface area across different lengths.

**Go/no-go gate — decide by 2026-10-13.** This stage has no slack under the pulled-in timeline,
and its duration is the one genuine unknown in the project (PANDORA's first run builds a local
MHC template database, on a disk-constrained machine). So the first task is a single peptide,
not the panel, and the decision is explicit rather than discovered late:

- Run PANDORA on **E01** alone and RMSD it against 1S9V, using the Checkpoint 2 apparatus and
  its measured scale. Record wall-clock time and disk used.
- **Accurate enough** (core RMSD within the ~1.16 Å cross-deposition ceiling) **and** ~50
  complexes fit the time and disk budget → proceed to the full panel.
- **Accurate but slow or disk-hungry** → scope down to a defensible subset: all 17 panel rows
  plus one matched control each (~34), dropping the second control tier. Say so in the report.
- **Not accurate enough** → that is a real limitation to state plainly, not a reason to reach
  for AlphaFold Server (see the constraint above). Fall back to interface-geometry scores that
  do not depend on sub-ångström accuracy, and lead the report with the sequence baseline as the
  primary comparison with structure as support.

- [x] PANDORA on E01; RMSD vs 1S9V; time + disk measured; gate decision recorded here.
      **2026-10-06: GATE PASSED — proceed to the full panel, with one protocol change.**

      Getting PANDORA to run at all took a compatibility layer: `csb-pandora` 0.9 is the only
      PyPI release and is broken two independent ways (three Biopython APIs removed since
      1.80; IMGT now https-only while PANDORA requests http, which kills both its structure
      download *and* its ID query). Full account in [`pandora_setup.md`](pandora_setup.md);
      `src/grains/pandora_compat.py` + `src/grains/pandora_db.py`, 30 tests.

      **The disk risk was unfounded.** PANDORA's own downloader would have blown past the
      free space (`gunzip` then `tar` leaves an uncompressed copy beside the extracted tree,
      from a 1.75 GB archive). Asking IMGT which PDB IDs are class II *first* and
      stream-extracting only those 205 members gives: **136 templates, 3.8 min, 0.24 GB**
      (`data/processed/pandora_db_build.json`). Modeling runs ~75 s–5 min per peptide for 20
      loop models, so ~40–50 complexes is one to four hours. **Neither time nor disk
      constrains this checkpoint.**

      **1S9V is itself one of the 136 templates**, carrying the exact peptide `LQPFPQPELPY`.
      So the gate could not be one run — `grains.validate_pandora` withholds templates in
      four conditions (`data/processed/pandora_gate.json`, `rmsd_pandora.csv`):

      | Condition | Template | State | Shared | Core RMSD |
      |---|---|---|---|---|
      | `full_db` (nothing withheld) | 1S9V | binary | 11/11 | **0.16–0.23 Å** ← circular |
      | `no_exact` (1S9V withheld) | 4OZI | ternary | 10/11 | 0.90–0.96 Å |
      | `no_core` (all E01-core templates withheld) | 6MFF | binary | 7/11 | **0.68–0.69 Å** |
      | `no_exact_binary` (1S9V + all ternary withheld) | 6MFG | binary | 10/11 | **0.56–0.58 Å** |

      Against the Checkpoint 2 scale: 0.24 Å within-crystal floor, 1.16 Å cross-deposition
      ceiling.

      *`full_db` quantifies the circularity.* 0.16–0.23 Å is **below the within-crystal
      floor** — closer to 1S9V than two copies in the same crystal are to each other. No
      prediction achieves that; it threaded the target onto the deposited answer. Reported so
      it can be dismissed, not quietly used.

      *The gate metric is `no_core`: 0.69 Å, well inside the 1.16 Å ceiling.* That is the
      condition resembling the real task — for a gluten epitope PANDORA usually has a close
      DQ2.5 template, but for a **control peptide from a safe grain it never will**, since no
      such complex has been crystallized. Accuracy holding up with only 7 of 11 residues
      shared is the result that licenses the full panel.

      **Unexpected finding, and it changes the protocol: template conformational state
      matters more than sequence similarity.** `no_core` (7/11 shared) *beat* `no_exact`
      (10/11) — 0.68 vs 0.90 Å. Cause: 4OZI is TCR-bound, and a homology model inherits its
      template's peptide conformation, so it reproduced 4OZI faithfully *including* 4OZI's own
      0.91–1.16 Å offset from binary 1S9V (measured independently in Checkpoint 2). The
      `no_exact_binary` condition isolates it with similarity held constant — 6MFG, binary,
      also 10/11 shared → **0.56–0.58 Å, a 0.35 Å improvement over the ternary template of
      equal similarity.**

      PANDORA ranks templates by peptide alignment score alone, which is blind to TCR
      presence, and **8 of the 12 DQ2.5 templates are ternary** — so left alone it will often
      pick one. This also contradicts the v1 scope: the project models peptide–MHC *without* a
      TCR, so a ternary template is the wrong reference state.

      **Protocol change, carried into the modeling step below:** model the panel against
      **binary templates only** (`modeling.load_database(binary_only=True)`, leaving 1S9V,
      6MFF, 6MFG, 6U3M). Caveat for the report: the effect is measured on **one peptide**, so
      it is consistent with Checkpoint 2's independent crystal comparison but it is n=1 — a
      panel peptide that models oddly is worth re-checking against the full template set.

      Also: record **the chosen template and its longest shared residue run per complex.**
      That is the covariate explaining each model's accuracy, and a toxic-vs-safe comparison
      in which toxic peptides systematically had better templates would be measuring template
      availability rather than binding.
- [ ] `build_controls.py` — length- and P/Q-matched control peptides from zein / kafirin /
      orzenin / avenin + a non-epitope gluten stretch → `data/processed/controls.csv`.
- [ ] Model toxic panel + controls + native pairs (~40–50 complexes) **with PANDORA**,
      **binary templates only** (gate finding above), recording the chosen template and
      its longest shared residue run for each complex.
- [ ] `src/grains/scoring.py` — per-complex: interface contact count, buried surface area
      (freesasa), PRODIGY ΔG / Kd, interface PAE + peptide pLDDT, anchor-pocket occupancy.
- [ ] `data/processed/scores.csv`; first comparison plot (toxic vs safe vs intermediate).

**Present:** score table + a boxplot/strip plot of one or two scores by celiac status.

---

## Checkpoint 4 — Analysis & baseline  (due 2026-10-30)

**Goal:** turn scores into a defensible answer.

- [ ] Sequence baseline: IEDB MHC-II API (`netmhciipan` method, alleles HLA-DQA1*05:01/
      DQB1*02:01 and DQ8) for every panel peptide → `data/processed/seq_baseline.csv`.
- [ ] `src/grains/deamidation.py` — TG2 QxP-type rules; check the `form` pairs are consistent.
- [ ] Stats: Mann–Whitney U (toxic vs safe) per score; effect sizes; deamidated vs native
      (paired). ROC / AUC using celiac status as the label.
- [ ] Compare structure scores vs the sequence baseline — where does structure add signal?

**Build a real null** (carried in from the ESM-2 comparison, 2026-10-06 — see
`background.md` §4). Two requirements, both cheap and both easy to skip by accident:

- **Permute the labels.** Shuffle toxic/safe a few thousand times, recompute the statistic, and
  report where the real value sits in that distribution. This is the honest version of "is the
  difference real" and it costs nothing once `scores.csv` exists.
- **Match decoys on length and P/Q composition, not just length.** Prolamins are
  compositionally extreme — mostly proline and glutamine — so generic decoy peptides are too
  easy to beat and would flatter the result. `build_controls.py` already matches on P/Q; the
  null has to as well.
- **Do not invent a p-value from a score.** `p = 1 - score` is not a p-value, is not uniform
  under the null, and any multiple-testing correction fed on it returns a confident-looking
  number with no guarantee behind it. Published instance of exactly this mistake noted in
  `background.md` §4.
- **Resample at the family level.** `E03`/`E06` are one observation (7 of 9 residues shared), so
  bootstrap over epitope families, not individual peptides, or the confidence intervals come out
  falsely narrow.

**Present:** 2–3 figures with statistics + a paragraph interpreting them.

---

## Checkpoint 5 — Complete draft  (due 2026-11-07)

**Goal:** the whole thing exists end to end. Rough in places is fine — complete is the bar, not
polished. Everything after this date improves a draft that already exists.

- [ ] `notebooks/analysis.ipynb` runs top to bottom from `data/processed/`.
- [ ] `report.md` (2–3 pp): question, methods, results, limitations, "what I'd do next".
      Introduction adapts `background.md` §1–2; the methods section's licensing paragraph
      adapts §3 and the Checkpoint 3 constraint.
- [ ] README updated with results summary + a headline figure.
- [ ] Repo tidy: docstrings, `pytest` green, `requirements.txt` complete, LICENSE.

**Present:** the repo + report as the portfolio piece.

---

## Finishing touches & fixes  (2026-11-07 → 11-28)

Three weeks against an existing draft. Nothing here is load-bearing — if a given item does not
happen, the portfolio piece is still complete and defensible. Ordered by value.

- [ ] Read the report cold after a few days away and fix what does not survive that. This is
      the single highest-value item in the window and the reason the draft moved to 11-07.
- [ ] Replace `results/figures/1S9V_dq25_groove.png` (Checkpoint 1's quick Biopython
      CA-backbone-trace render) with a real ChimeraX or PyMOL cartoon: helices as ribbons,
      beta-sheet as arrows, peptide as sticks colored by atom, viewed straight down into the
      groove. The current figure only shows the peptide is roughly in the right place; it
      can't show anchor-pocket chemistry (the deamidated Glu side chains dropping into the
      DQ2.5 pockets), which is the actual story worth a headline figure. Needs disk space
      freed up first (no Xcode CLT / low disk currently blocks a ChimeraX/PyMOL install).
- [ ] Application text itself — the research-interest paragraphs that reference this work.
- [ ] Open Checkpoint 1 item: the `E10` `DQ8-glia-γ1a` vs `g1` name-suffix question, if the
      Sollid 2020 text becomes reachable (paywalled as of 2026-10-05).
- [ ] Optional stretch: a second HLA allele (DQ8), or ancient-wheat (einkorn/emmer) prolamins.
- [ ] Optional stretch: run the oat avenin-E candidates from the ESM-2 study
      (`background.md` §4) through PANDORA + PRODIGY and place them against the panel's
      toxic/safe score distributions — an independent structural read on candidates a sequence
      classifier cannot adjudicate. Needs the cultivar sequences; that repo is private.

---

## Risk log

| Risk | Mitigation |
|------|-----------|
| AlphaFold Server rate limits / queue | PANDORA (CPU) is the fallback modeling route; start submissions early |
| Sparse UniProt coverage for teff / millets | Drop to a smaller grain set; note the data gap explicitly in the report |
| Peptide sequence errors propagate | Checkpoint 1 verification gate; `verified == yes` required before modeling |
| Scope creep | TCR step, MD refinement, and extra alleles are explicitly out of scope for v1 |
| **Checkpoint 3 has no slack** (pulled-in timeline) | Explicit go/no-go gate at 2026-10-13 on a single peptide, with a pre-agreed scoped-down panel (~34 complexes) and a pre-agreed fallback if accuracy is poor. The decision is made on measured time and disk, not discovered at 10-20 |
| **PANDORA template database fills the disk** | Measure disk used during the E01 gate run before committing to the panel; machine was at 94% / 13 GB free on 2026-10-05 |
| Compressed Checkpoint 4/5 squeezes the writeup | The 11-07 → 11-28 window is revision against a complete draft, so the squeeze costs polish, not completeness |
