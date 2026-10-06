# CLAUDE.md

Orientation for a session picking this repo up. Two linked projects live here.

## 1. Research project (repo root) — grain peptides × HLA-DQ2.5

**Question:** do grain-derived peptides that bind HLA-DQ2.5 explain which cereals are toxic in
celiac disease? Model known gluten epitopes + matched control peptides from safe grains in the
DQ2.5 groove, score binding, test whether it separates toxic (wheat/barley/rye) from safe.
This is the user's comp-bio PhD-application portfolio piece (applications due end of Nov 2026).

**Framing that keeps getting muddled — hold it:**
- No "gluten receptor". Target = the **HLA-DQ2.5 peptide-binding groove** (MHC class II).
- Model short **peptides** (~9–17 aa), not whole proteins.
- Genome models (Evo2 etc.) are the wrong tool — this is peptide–MHC immunology. Use
  AlphaFold Server / PANDORA for structure, PRODIGY for scoring, IEDB API for a sequence
  baseline. All CPU-only friendly.

**Layout**
```
data/reference/    curated, version-controlled inputs
  grains.csv          16 grains: species, NCBI taxon id, storage-protein family, toxic/safe
  peptides.csv        17 gluten T-cell epitopes; `verified` column tracks checking status
  peptides_NOTES.md   verification checklist against Sollid et al. 2012 / 2020 + IEDB
  pdb_structures.md   candidate HLA-DQ2.5 crystal structures (e.g. 4OZF) to validate against
data/raw/          downloaded proteomes (git-ignored; regenerate with proteomes.py)
data/processed/    proteome_index.csv, proteome_summary.csv
src/grains/        config.py, proteomes.py, fetch_hla.py, render_structure.py,
                   superpose.py (RMSD), modeling.py (AlphaFold/PANDORA routes)
tests/             pytest smoke tests for the reference data + structure handling
docs/background.md    the biology, the question, what each pipeline stage is for — written to
                      be read cold; becomes the report's introduction
docs/checkpoints.md   the 5-checkpoint plan, dated 2026-09-08 → 2026-11-28
```

**Explaining this project to the user:** they are learning the immunology alongside the code and
have asked twice for a simpler, deeper explanation — pitch at the level of `docs/background.md`,
which was written for exactly that and is the shared reference. Define terms on first use, build
from the biology up rather than from the tooling down, and don't assume MHC/peptide vocabulary.

**Run**
```bash
pip install -r requirements.txt          # adds biopython to the Anaconda base env
python -m grains.proteomes               # fetch seed-storage proteins from UniProt
python -m pytest -q                      # 50 tests, currently green
```
Python is Anaconda at `/opt/anaconda3/bin/python3` (3.13). `PYTHONPATH=src` for the module,
or rely on `pyproject.toml`'s `pythonpath`.

**State (2026-10-06):** **Checkpoints 1 and 2 both done.** Checkpoint 3 (model + score the
panel) is next, due 2026-10-20.

**Timeline pulled in 2026-10-06 at the user's request:** everything substantive done by
**2026-11-07**, then 11-07 → 11-28 is finishing touches and fixes. CP3 keeps 10-20; CP4 moved
11-03 → **10-30**; CP5 moved 11-17 → **11-07** and rescoped to *complete draft* rather than
polished. The squeeze lands on CP4/CP5 and buys three weeks of revision against an existing
draft instead of eleven days of slippage buffer. Consequence to respect: **CP3 now has no
slack**, so it carries an explicit go/no-go gate at **2026-10-13** — PANDORA on E01 alone,
measuring accuracy *and* wall-clock *and* disk, with a pre-agreed scoped-down panel (~34
complexes) and a pre-agreed fallback if accuracy is poor. Don't let CP3 drift past that gate.

*Checkpoint 1.* Reference data built; proteome fetcher run (1,595 proteins; coverage gap for
teff/pearl-millet/amaranth is real, note it, don't hide it). All 17 peptide rows are now
verified, IEDB-identified (`iedb_id`, 12/15 cored rows) and mapped to a UniProt
accession+range — and **gene-faithful**: 13 deamidated rows each need ≥1 Q→E against the gene,
4 native rows match literally. Four classes of silent error were found and fixed along the
way: five wrong PMIDs across two audits (incl. `22371225`, cited by *every* epitope row);
**E09** carrying a native `extended_sequence` under a `deamidated` label; and **five rows with
fabricated flanks** whose cores were fine but whose surrounding residues matched no gene.
Only open item: the E10 `γ1a`/`g1` naming question, which needs the paywalled Sollid 2020.

*Checkpoint 2 — the pipeline reproduces a known complex.* Validation structure is **`1S9V`**,
not 4OZF (all of 4OZF–4OZI are TCR-bound); its peptide is panel epitope E01 exactly. The
AlphaFold Server model hits **core RMSD 0.59–0.73 Å**, better than two independent crystals of
the same complex agree with each other (0.91–1.16 Å), against a 0.24 Å within-crystal floor.
Checked for circularity rather than assumed: 1S9V wasn't a template, 4OZI was, and the model
is equidistant from both rather than hugging its template. Figure in
`results/figures/af_vs_1S9V_overlay.png`. **MODELLER + PANDORA are installed and licensed.**

*Standing constraint:* Checkpoint 3 locks **PANDORA as the only route** for complexes that get
PRODIGY-scored (AlphaFold Server ToS bars its output from binding predictors);
`modeling.assert_scorable()` enforces it in code.

**The recurring failure mode in this project is the silently plausible wrong number** — a
result that arrives looking computed rather than erroneous. Paid for so far: a disordered
7-residue loop gap in 1S9V chain B that a naive CA walk splices across; **4OZI's peptide chains
being the reverse of their chain lettering** (J binds A/B — pairing by letter gives a 62 Å
peptide RMSD); wrong PMIDs that still resolve to real papers; a readiness check that passed on
an unvalidated MODELLER license. The habit that catches these: assert the assumption in code
and keep the bad case as a regression test. 50 tests, green.

Two carry-forwards into Checkpoint 3/4 statistics: **E03 and E06 overlap by 7 of 9 residues**
(wheat ω-gliadin vs rye secalin — the same motif, which is itself part of why rye is toxic) so
they are not independent observations; and panel peptide lengths run 12–18 aa, so length is a
confound for any raw interface score.

**`git` is initialized**
(2026-09-14, initial commit `d3e9ba8`) — this Mac has no Xcode Command Line Tools and is
low on disk space, so git commands use GitHub Desktop's bundled binary instead:
`/Applications/GitHub Desktop.app/Contents/Resources/app/git/bin/git`. Alias or invoke by
full path; the system `git` at `/usr/bin/git` still just triggers the CLT install prompt.
That bundled git is a relocated build: invoked directly it resolves its own support dirs to
`//libexec/git-core` and `//share/git-core/templates`, so it can't find its helper
subcommands, falls back to `PATH`, hits the `/usr/bin/git` stub, and pops the CLT installer
mid-commit (the commit still succeeds). `GIT_EXEC_PATH` and `GIT_TEMPLATE_DIR` are now set
to the bundle's real `libexec/git-core` and `share/git-core/templates` in
`~/.claude/settings.json` `env`, which silences it. Outside Claude Code — a bare terminal,
a script — set them there too, or the prompt comes back.

## 2. Teaching project — `teaching/data8-project/`

A five-part intro-data-science project for the user's Sierra College Tutor Center role:
*can you spot a gluten grain by counting amino acids?* Scoped-down slice of project 1.

- `INSTRUCTOR.md` — full spec: 5 parts, per-question answers, review guidance, common mistakes.
- `TUTOR_PROJECT.md` — rubric-facing rationale (audience, purpose, CRLA alignment, APA sources).
- `datascience_vs_pandas.md` — side-by-side reference for the two notebook tracks.
- `build_dataset.py` → `data/grains.csv` (251 proteins) + `data/grains_mystery*.csv` (18).
- `build_notebooks.py` → the 4 `.ipynb` (datascience/pandas × project/solutions).

**Conventions:** prediction-focused (Part 5). Two parallel tracks (`datascience` +
`pandas`) generated from one source file. Judged on the held-out **mystery set**, not a
random split (random split is near-saturated). A `check_*` helper after every question,
including written ones.

**State (2026-10-05):** spec + dataset + rationale done. Parts 0–3 built in both tracks;
solutions execute clean. Parts 0–2 reviewed and revised; **Parts 4–5 not built yet.**
(Bootstrap part cut; parts run 0–5. The module is unscored.) Regenerate *and verify* with
`./rebuild.sh`, which runs `build_notebooks.py` and then executes both solutions notebooks —
Part 3's simulation makes the `datascience` run take ~11 s, which is expected.

## Conventions (both projects)

- **American spelling everywhere** — code, comments, docstrings, prose.
- Teaching code: explicit step-by-step, not clever one-liners.
- The user reviews closely and iterates; build in visible stages, explain terms on first use.
- Attribution for commits/PRs is set in the session; end commit messages with the
  `Co-Authored-By` / `Claude-Session` trailer already configured.

## Published

Visual overview of both projects (private artifact, user shares manually):
https://claude.ai/code/artifact/0af628d9-1263-452f-9698-e21c061b44b8
