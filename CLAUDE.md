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
docs/checkpoints.md   the 5-checkpoint plan, dated 2026-09-08 → 2026-11-28
```

**Run**
```bash
pip install -r requirements.txt          # adds biopython to the Anaconda base env
python -m grains.proteomes               # fetch seed-storage proteins from UniProt
python -m pytest -q                      # 27 tests, currently green
```
Python is Anaconda at `/opt/anaconda3/bin/python3` (3.13). `PYTHONPATH=src` for the module,
or rely on `pyproject.toml`'s `pythonpath`.

**State (2026-10-05):** Checkpoint 1 done, Checkpoint 2 mostly done (due 10-06). Reference
data built; proteome fetcher run (1,595 proteins; coverage gap for rye/teff/pearl-millet/
amaranth is real, note it, don't hide it). All 17 peptide rows verified against Sollid 2012
(two wrong PMIDs caught and fixed); still open are the `iedb_id` and UniProt accession+range
columns, the Sollid 2020 update, and a direct Tye-Din 2010 check for N01–N03 (verified by
derivation only). Validation structure is **`1S9V`**, not 4OZF (all of 4OZF–4OZI are
TCR-bound); its peptide is panel epitope E01 exactly.

Checkpoint 2: `fetch_hla.py` extracts the DQ2.5 ectodomains from 1S9V's SEQRES records;
`superpose.py` measures the **crystal-vs-crystal RMSD baseline** — E01 core 0.24 Å between
copies in one crystal, 1.16 Å between independent depositions, so a model within ~1.2 Å has
matched experimental reproducibility. Two things are blocked on manual steps: the AlphaFold
Server upload (no public API; submission JSON is written and ready) and the **MODELLER
license key**, which blocks PANDORA — the critical path, since Checkpoint 3 locks PANDORA in
as the only route for complexes that get PRODIGY-scored. `modeling.assert_scorable()` enforces
that in code.

Two traps already paid for, both of the "silently plausible wrong number" kind — assert the
assumption rather than trusting the file: a disordered 7-residue loop gap in 1S9V chain B that
a naive CA walk splices across, and **4OZI's peptide chains being the reverse of their chain
lettering** (J binds A/B, I binds C/D; pairing by letter gives a 62 Å peptide RMSD).

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
