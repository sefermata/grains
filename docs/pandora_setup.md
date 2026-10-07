# Getting PANDORA to run (2026-10-06)

PANDORA is the modeling route for every complex that gets PRODIGY-scored, so Checkpoint 3
cannot start without it. Installing it took more than `pip install`: **`csb-pandora` 0.9 is the
only release on PyPI, it was last built against a 2022-era environment, and both its Biopython
calls and its data downloads are now broken.** None of it is unfixable, but none of it is
guessable either, so this records what was wrong and what was done.

Everything here is encoded in `src/grains/pandora_compat.py` (the shims) and
`src/grains/pandora_db.py` (the data acquisition), with tests in `tests/test_pandora.py`.
Nothing in `site-packages` was edited — all of it is patched at import time, so a reinstall of
PANDORA does not silently undo the fixes.

---

## 1. Install

```bash
conda install -c salilab modeller            # needs a free academic license key
pip install csb-pandora                      # PANDORA 0.9 (internal version string: 1.0)
conda install -c bioconda muscle=3.8.1551    # NOT muscle 5.x -- see §4
```

MODELLER's license key goes in its `config.py`. Note that `config.py` takes **precedence over
the `KEY_MODELLER` environment variable**, so leaving the placeholder `XXXX` in the file defeats
a correctly set env var. `grains.modeling.pandora_available()` validates the license by actually
constructing a `modeller.Environ()` rather than merely importing the module — an import succeeds
on an unlicensed install, which is how an earlier readiness check reported "ready" wrongly.

---

## 2. Three removed Biopython APIs

PANDORA was written against Biopython < 1.80 and calls three things that no longer exist. This
project runs Biopython 1.88, and the import fails at module load, so *nothing* ran — not even
the database build.

| PANDORA calls | Removed in | Reachable from |
|---|---|---|
| `Bio.PDB.Polypeptide.three_to_one` | 1.80 | `PMHC/Model.py` `renumber()` |
| `Bio.pairwise2` | 1.84 | `renumber()`, NetMHCpan anchor prediction |
| `Bio.Align.Applications.MuscleCommandline` | 1.82 | `Pandora/Align.py`, every MSA call |

**Downgrading Biopython was rejected.** The Checkpoint 2 validation — the measured RMSD scale
that makes every later number interpretable — was produced on 1.88, and 50 tests were green
against it. Trading a verified result for an unverified dependency to satisfy an unmaintained
package is the wrong direction.

Instead `pandora_compat` injects the three names. The important design decision is that **two
of the three raise `NotImplementedError` instead of being reimplemented:**

`three_to_one` and `pairwise2` are reachable only from `Model.renumber()` (a benchmarking
helper) and from NetMHCpan-based anchor prediction. This project passes `benchmark=False` and
supplies anchors explicitly, so neither runs. Reimplementing them would ship an untested
approximation of an alignment into a structural pipeline on a path nothing exercises — which is
exactly the *silently plausible wrong number* this project keeps paying for. A loud crash on an
unexpected path is strictly better. Both stubs name their correct replacement in the error
message, so if a future change does reach one, the fix is a two-line edit plus a test.

`MuscleCommandline` is different: it genuinely runs (see §4), so it genuinely works.

---

## 3. IMGT is https-only, and PANDORA asks for http

Independent of Biopython, PANDORA's whole data-acquisition path is dead. `www.imgt.org` no
longer serves plain HTTP — and not as a redirect, it returns nothing at all (`curl` exit code
000), so there is nothing for `urllib` to follow.

| PANDORA calls | Problem |
|---|---|
| `download_unzip_imgt_structures` | `http://` URL, and shells out to `wget`, absent on macOS |
| `download_ids_imgt` | `http://` URL — **and this one runs even with `download=False`** |

`download_ids_imgt` is replaced with an https equivalent in `pandora_compat`, parsing the
response exactly as the original did so the ID list is byte-identical to what PANDORA would have
produced when its URL worked. The patch is applied to the module object and then **asserted** to
have taken effect, because `Database.py` resolves the function off `Database_functions` at call
time — a future PANDORA using `from ... import download_ids_imgt` would bind the old function
directly and quietly keep using http.

### The disk problem, and why it evaporated

PANDORA's downloader does `gunzip archive.tgz` then `tar -xvf archive.tar`, which leaves a fully
uncompressed tar *beside* the extracted tree. The archive is **1.75 GB compressed**; that peak
did not fit the free space on this machine, and it was the reason Checkpoint 3 was flagged as a
disk risk at all.

`grains.pandora_db` avoids it with two reductions:

- **Class II only.** This project models HLA-DQ2.5. Parsing class I too would roughly double
  time and disk for templates that can never be selected.
- **Extract only the class II members.** IMGT's flat-file set covers every immune structure it
  holds — antibodies, TCRs, MHC class I — as one `IMGT-<id>.pdb.gz` per entry. Asking IMGT
  *first* which PDB IDs are class II (205 of them) means extracting 205 members rather than
  thousands. The archive is streamed with `tarfile.open(..., "r|gz")`, so the uncompressed tar
  is never written at all, and the archive is deleted afterward.

Measured result (`data/processed/pandora_db_build.json`):

| | |
|---|---|
| IMGT class II structures listed | 205 |
| Members extracted | 203 (6CQN and 6CQQ are listed but absent from the archive) |
| Templates successfully parsed | **136** |
| Wall clock | **3.8 min** |
| Disk used | **0.24 GB** |

**So the disk risk was unfounded once the acquisition was done properly** — a quarter of a
gigabyte, not tens. The 67 structures that fail to parse are PANDORA's own quality filtering
(non-canonical residues in the peptide chain, and so on), not an error here.

PANDORA also assumes a directory tree it never ships, and `os.chdir`s into those paths without
creating them, so a missing directory surfaces as a bare `FileNotFoundError`. `ensure_tree()`
creates all seven, including `outputs/` — which is needed only at the very *end* of a
successful model run, where `Pandora.__log` writes its logfile. That one crashed a model that
had already completed successfully.

---

## 4. MUSCLE 3, specifically

`Align.py` aligns the MHC class II M and N chains against the templates with MUSCLE, and uses it
as:

```python
msl = MuscleCommandline(input=..., out=...)
os.system(str(msl) + ' -quiet')
```

So the shim only has to render a runnable command line. It emits Biopython's **MUSCLE v3 syntax**
(`muscle -in IN -out OUT`), which is why the pinned version matters: **MUSCLE 5 renamed those
flags to `-align`/`-output`**. Under v5 the command would exit without writing the output file,
and PANDORA — which does not check — would carry on with a missing alignment. `muscle_syntax_ok()`
checks the installed version's banner, and a test asserts the rendered command contains `-in`/
`-out` and *not* `-align`/`-output`.

---

## 5. The allele name is not guessable

PANDORA keys class II alleles with IMGT's spelling as it stores it:

```
HLA-DQA1*0501      correct   -- asterisk kept, four digits, no colon
HLA-DQA1*05:01     wrong     -- the standard nomenclature
HLA-DQA10501       wrong     -- the first guess made here
```

This matters more than a typo normally would, because **PANDORA does not reject an unknown
allele name.** It scores every template and returns the best one, so a misspelling produces a
model built on some other allele's structure — a plausible-looking result with no error
anywhere. `modeling.check_allele_names()` validates the names against the built database and
reports the closest real ones; it caught the wrong first guess immediately, which is the whole
argument for having it in the call path rather than trusting the constant.

Of the 136 templates, 12 carry a DQ2.5 chain.

---

## 6. Anchors are derived, never hardcoded

MHC class II binds a 9-residue core through four anchor pockets, **P1/P4/P6/P9** — core offsets
0, 3, 5, 8. PANDORA wants them as 1-based positions in the *full* peptide, so
`modeling.class_ii_anchors()` adds back the core's offset within the extended sequence.

It is computed per row rather than fixed, because the panel's cores sit at different offsets:
E01's core starts at position 3 of `LQPFPQPELPY`, E04's at position 1 of `PQPEQPFPWQPEQ`. A
single hardcoded anchor list would be right for one row and quietly wrong for the next. A test
asserts the offsets actually differ across the panel, so the derivation cannot become moot
without failing.

For E01 this gives anchors `[3, 6, 8, 11]` → residues `P, Q, E, Y`, which puts the
TG2-derived glutamate on **P6** — one of DQ2.5's positively charged pockets, exactly as the
biology predicts (see [`background.md`](background.md) §1, Link 6). That agreement is a check on
the anchor convention, not a coincidence.

One subtlety worth keeping straight: **P7 is a charge-preferring pocket for DQ2.5 but is not one
of PANDORA's four geometric anchors.** Panel row E08 (`PQQSFPEQQ`) is where that bites — its
only glutamate is at P7. The test for "every deamidated DQ2.5 core carries a glutamate at a
charge pocket" therefore checks P4/P6/**P7**, not the anchor set.

---

## 7. 1S9V is in the template database

Found while watching the build log: `Parsing 1S9V`. The Checkpoint 2 validation structure is one
of PANDORA's 136 templates, carrying the exact peptide `LQPFPQPELPY`.

So asking PANDORA to model the validation peptide against the full database hands it the answer:
it selects 1S9V and threads the target onto it. The RMSD would be excellent and would mean
nothing — the same circularity Checkpoint 2 had to reason about for AlphaFold, except here the
template set is ours to control, so it can be settled by withholding rather than inferred.

`grains.validate_pandora` therefore runs the gate under three conditions:

| Condition | Withheld | Measures |
|---|---|---|
| `full_db` | nothing | the circular ceiling, reported so it can be dismissed |
| `no_exact` | 1S9V | can PANDORA rebuild the complex without the complex? |
| `no_core` | 1S9V, 4OZI, 6MFG (every template whose peptide contains `PFPQPELPY`) | **the honest estimate** |

`no_core` is what the gate decision rests on, because it resembles the real task. For a gluten
epitope PANDORA will usually have a close DQ2.5 template; for a **control peptide from a safe
grain it will not** — no such complex has ever been crystallized. If accuracy collapses once
near-identical templates are withheld, the toxic-vs-safe comparison would partly be measuring
template availability rather than binding, which is a confound severe enough to matter to the
design.

Each run also records the chosen template and the longest contiguous residue run it shares with
the target, since that is the covariate that explains the RMSD: a template sharing 10 of 11
residues is nearly the answer, one sharing 4 is a real prediction.

A fourth condition, `no_exact_binary`, withholds 1S9V **and every ternary template**, which
leaves 6MFG — binary, and sharing 10 of 11 residues exactly as 4OZI does. Similarity held
constant, conformational state varied, which turns an observation into a controlled comparison.

### Result

| Condition | Template | State | Shared | Core RMSD vs 1S9V |
|---|---|---|---|---|
| `full_db` | 1S9V | binary | 11/11 | **0.16–0.23 Å** (circular) |
| `no_exact` | 4OZI | ternary | 10/11 | 0.90–0.96 Å |
| `no_core` | 6MFF | binary | 7/11 | **0.68–0.69 Å** |
| `no_exact_binary` | 6MFG | binary | 10/11 | **0.56–0.58 Å** |

Floor 0.24 Å, ceiling 1.16 Å. **Gate metric is `no_core` → 0.69 Å → PASS.**

`full_db` lands *below the within-crystal floor*, which is the signature of copying rather than
predicting — exactly why it is reported instead of used.

The unexpected part: **`no_core` (7/11 shared) beat `no_exact` (10/11)**. Template conformational
state matters more than sequence similarity here. 4OZI is TCR-bound, and a homology model
inherits its template's peptide conformation, so building on it reproduced 4OZI faithfully
*including* 4OZI's own 0.91–1.16 Å offset from binary 1S9V — the offset Checkpoint 2 had already
measured from the crystals alone. `no_exact_binary` confirms it: at equal similarity, the binary
template is **0.35 Å better**.

Since PANDORA ranks templates by peptide alignment score alone — blind to TCR presence — and 8
of the 12 DQ2.5 templates are ternary, Checkpoint 3 models the panel against **binary templates
only** (`modeling.load_database(binary_only=True)`; see `modeling.binary_templates_only`). That
also agrees with the v1 scope, which deliberately excludes the TCR.

Measured on one peptide, so n=1 — consistent with Checkpoint 2's independent crystal comparison,
but a panel peptide that models oddly is worth re-checking against the full template set.

Results: `data/processed/rmsd_pandora.csv` and `data/processed/pandora_gate.json`.
