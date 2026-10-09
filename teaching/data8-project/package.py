"""Build the delivery bundles for the Tutor Center.

Produces two zips in dist/:

  gluten-grains-module-<date>.zip    what you hand the Tutor Center. Everything:
                                     both student notebooks, both solutions (executed,
                                     so the outputs are readable without running
                                     anything), the data, the regeneration source, and
                                     a copy of the student zip ready to forward.

  gluten-grains-students-<date>.zip  what students get. The two student notebooks and
                                     the data they need, and nothing else. No
                                     solutions, and no mystery-set answer key.

INSTRUCTOR.md and TUTOR_PROJECT.md are deliberately left out of both.

Usage:  python package.py          (add --skip-execute to build without running the
                                    solutions, which is much faster but ships them blank)
"""

from __future__ import annotations

import argparse
import datetime as dt
import shutil
import subprocess
import sys
import zipfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
DIST = HERE / "dist"
STAGE = DIST / "_staging"

STUDENT_NOTEBOOKS = ["project_datascience.ipynb", "project_pandas.ipynb"]
SOLUTION_NOTEBOOKS = ["solutions_datascience.ipynb", "solutions_pandas.ipynb"]
# Generated, and needed at runtime beside the notebooks - not source/ material.
CHECKS_MODULES = ["checks_datascience.py", "checks_pandas.py"]
STUDENT_DATA = ["grains.csv", "grains_mystery.csv"]
INSTRUCTOR_ONLY_DATA = ["grains_mystery_key.csv"]
SOURCE_FILES = ["build_dataset.py", "build_notebooks.py", "rebuild.sh", "package.py"]

# Anything matching these must never reach the student zip.
FORBIDDEN_IN_STUDENT = ["solutions_", "grains_mystery_key", "INSTRUCTOR", "TUTOR_PROJECT"]

REQUIREMENTS = """\
# Install with:  pip install -r requirements.txt
datascience>=0.17
pandas>=2.0
numpy>=1.24
matplotlib>=3.7
notebook>=7.0
"""


def run(cmd: list[str], **kw) -> subprocess.CompletedProcess:
    return subprocess.run(cmd, cwd=HERE, capture_output=True, text=True, **kw)


def execute_solutions(into: Path) -> list[str]:
    """Run each solutions notebook and write the executed copy into `into`.

    Shipping solutions with their outputs means a tutor can read the expected
    numbers and plots without installing anything first.
    """
    problems = []
    for name in SOLUTION_NOTEBOOKS:
        print(f"    executing {name} ...", end=" ", flush=True)
        proc = run([
            "jupyter", "nbconvert", "--to", "notebook", "--execute", name,
            "--output", str(into / name),
        ])
        if proc.returncode == 0:
            print("ok")
        else:
            print("FAILED")
            problems.append(f"{name}: {proc.stderr.strip().splitlines()[-1] if proc.stderr.strip() else 'unknown error'}")
            shutil.copy2(HERE / name, into / name)
    return problems


def build_student_tree(root: Path) -> None:
    root.mkdir(parents=True, exist_ok=True)
    for name in STUDENT_NOTEBOOKS + CHECKS_MODULES:
        shutil.copy2(HERE / name, root / name)
    shutil.copy2(HERE / "datascience_vs_pandas.md", root / "datascience_vs_pandas.md")
    data = root / "data"
    data.mkdir(exist_ok=True)
    for name in STUDENT_DATA:
        shutil.copy2(HERE / "data" / name, data / name)
    (root / "requirements.txt").write_text(REQUIREMENTS)
    (root / "README.md").write_text(STUDENT_README)


def build_module_tree(root: Path, skip_execute: bool) -> list[str]:
    root.mkdir(parents=True, exist_ok=True)
    for name in STUDENT_NOTEBOOKS + CHECKS_MODULES:
        shutil.copy2(HERE / name, root / name)
    shutil.copy2(HERE / "datascience_vs_pandas.md", root / "datascience_vs_pandas.md")

    problems = []
    if skip_execute:
        for name in SOLUTION_NOTEBOOKS:
            shutil.copy2(HERE / name, root / name)
    else:
        problems = execute_solutions(root)

    data = root / "data"
    data.mkdir(exist_ok=True)
    for name in STUDENT_DATA + INSTRUCTOR_ONLY_DATA:
        shutil.copy2(HERE / "data" / name, data / name)

    source = root / "source"
    source.mkdir(exist_ok=True)
    for name in SOURCE_FILES:
        shutil.copy2(HERE / name, source / name)
    (source / "README.md").write_text(SOURCE_README)

    (root / "requirements.txt").write_text(REQUIREMENTS)
    (root / "README.md").write_text(MODULE_README)
    return problems


def zip_tree(tree: Path, target: Path, arc_root: str) -> None:
    with zipfile.ZipFile(target, "w", zipfile.ZIP_DEFLATED) as z:
        for path in sorted(tree.rglob("*")):
            if path.is_file():
                z.write(path, Path(arc_root) / path.relative_to(tree))


def verify_student_zip(path: Path) -> list[str]:
    """The answer key must not be in the student bundle. Check, don't assume."""
    names = zipfile.ZipFile(path).namelist()
    return [n for n in names if any(bad in n for bad in FORBIDDEN_IN_STUDENT)]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--skip-execute", action="store_true",
                    help="don't run the solutions notebooks (faster, but they ship with no outputs)")
    args = ap.parse_args()

    stamp = dt.date.today().isoformat()
    if DIST.exists():
        shutil.rmtree(DIST)
    STAGE.mkdir(parents=True)

    print("\n==> Regenerating notebooks from build_notebooks.py")
    proc = run([sys.executable, "build_notebooks.py"])
    if proc.returncode != 0:
        print(proc.stdout + proc.stderr)
        print("!!! build_notebooks.py failed - nothing was packaged.")
        return 1
    print("    " + proc.stdout.strip().replace("\n", "\n    "))

    print("\n==> Building the student bundle")
    student_tree = STAGE / "gluten-grains"
    build_student_tree(student_tree)
    student_zip = DIST / f"gluten-grains-students-{stamp}.zip"
    zip_tree(student_tree, student_zip, "gluten-grains")
    leaks = verify_student_zip(student_zip)
    if leaks:
        print("!!! the student zip contains files it must not:")
        for n in leaks:
            print("      " + n)
        return 1
    print(f"    {student_zip.name}  (no solutions, no answer key - verified)")

    print("\n==> Building the Tutor Center bundle")
    module_tree = STAGE / "gluten-grains-module"
    problems = build_module_tree(module_tree, args.skip_execute)

    handout = module_tree / "handout"
    handout.mkdir(exist_ok=True)
    shutil.copy2(student_zip, handout / student_zip.name)

    module_zip = DIST / f"gluten-grains-module-{stamp}.zip"
    zip_tree(module_tree, module_zip, "gluten-grains-module")
    print(f"    {module_zip.name}")

    shutil.rmtree(STAGE)

    print()
    if problems:
        print("!!! solutions did not execute cleanly - they shipped without outputs:")
        for p in problems:
            print("      " + p)
    size = lambda p: f"{p.stat().st_size / 1024:.0f} KB"
    print(f"Done. In {DIST.relative_to(HERE.parent.parent)}/")
    print(f"  {module_zip.name:45} {size(module_zip):>8}   <- send this one")
    print(f"  {student_zip.name:45} {size(student_zip):>8}   <- also inside it, under handout/")
    print()
    return 1 if problems else 0


MODULE_README = """\
# Spotting Gluten Grains from Protein Sequence

A self-contained data-science project for Data 8 students: *can the amino-acid makeup of a
protein tell a gluten grain from a safe one?* Students build features from raw protein
sequence, test whether the difference they find is real, engineer a sharper feature, and
end with a classifier they judge on proteins it has never seen.

No biology background is needed, by students or by whoever runs it. Everything required is
stated in Part 0 of the notebook.

## What to hand out

`handout/gluten-grains-students.zip` is the student version: the two notebooks and the data,
with **no solutions and no answer key**. That is the file to forward or upload.

Keep this bundle for yourself. It additionally contains the worked solutions and
`data/grains_mystery_key.csv`, which holds the labels for the held-out mystery set.

## What is here

| | |
|---|---|
| `project_datascience.ipynb` | student notebook, Berkeley `datascience` library |
| `project_pandas.ipynb` | student notebook, `pandas` |
| `solutions_datascience.ipynb` | worked solutions, saved with their outputs |
| `solutions_pandas.ipynb` | worked solutions, saved with their outputs |
| `checks_*.py` | the `check_` helpers the notebooks import. Generated; holds expected answers |
| `datascience_vs_pandas.md` | side-by-side reference for every operation the project uses |
| `data/` | the dataset (see below) |
| `source/` | scripts that regenerate everything, if the module ever needs editing |

The solutions are saved **with their outputs**, so you can read the expected numbers and
plots without installing or running anything.

## The two tracks

The same project, twice. `project_datascience.ipynb` uses Berkeley's `datascience` library,
the one Data 8 teaches with. `project_pandas.ipynb` uses `pandas`, the one used outside the
classroom. The questions, the data, and the answers are identical; only the library changes.
Students pick one. `datascience_vs_pandas.md` maps each operation to the other, and two
short exercises ask students to translate between them.

## Running it

```bash
pip install -r requirements.txt
jupyter notebook
```

Then open whichever notebook you want. Run the cells in order from the top — later cells
depend on variables set by earlier ones.

**Keep the `data/` folder and the `checks_*.py` files next to the notebooks.** The notebooks
load `data/grains.csv` by relative path and import their checks from the matching module, so
moving a notebook out of this folder on its own will break it.

The checks were pulled out of the notebook deliberately. Inline they ran to 260 lines a
student had to scroll past before the first question, and the hints quote expected answers.
Out here they still work exactly the same and are no longer in the way.

If local installation is not an option, the notebooks will run in Google Colab, but you have
to upload the `data/` folder alongside the notebook and run `!pip install datascience` in the
first cell. Nothing else changes.

## Self-checks

Every question ends in a `check_` call that prints `OK` or a specific hint. They are
guidance, not grading — there is no score anywhere in this module, and nothing is recorded.
The written questions are checked too, but only for "you actually wrote something".

## The data

`data/grains.csv` — 251 seed-storage proteins from 10 grains, 118 from gluten grains
(wheat, barley, rye) and 133 from gluten-free ones. One row is one protein. The only real
input is the `sequence` column; everything else students use, they compute.

`data/grains_mystery.csv` — 18 unlabeled proteins students classify at the end.
Deliberately atypical ones, which is the point: on a random split of the main dataset a
composition-based classifier is already near-perfect, so it cannot tell a good feature set
from a careless one. The mystery set can.

`data/grains_mystery_key.csv` — **the answers to the mystery set. Instructor copy only.**
It is not in the student zip.
"""

STUDENT_README = """\
# Spotting Gluten Grains from Protein Sequence

Wheat, barley, and rye trigger celiac disease. Rice, maize, oats, quinoa, millet, and
buckwheat do not. This project asks whether you can tell which is which from the protein
sequence alone — by counting amino acids.

**No biology background needed.** Everything you need is explained in Part 0.

## Pick one notebook

| | |
|---|---|
| `project_datascience.ipynb` | uses the Berkeley `datascience` library, the one from class |
| `project_pandas.ipynb` | uses `pandas`, the one used outside class |

Same project, same questions, same answers — only the library differs. Start with the
`datascience` one if you took Data 8. `datascience_vs_pandas.md` shows every operation
written both ways, side by side.

## Running it

```bash
pip install -r requirements.txt
jupyter notebook
```

Then open your notebook and run the cells in order from the top. Later cells use variables
that earlier ones create, so skipping around will give you errors that are not your fault.

**Keep the `data/` folder and the `checks_*.py` files next to the notebook.** The notebook
loads `data/grains.csv` by relative path, and the setup cell imports its checks from the
matching `checks_` file. Move the notebook somewhere on its own and neither will be found.

## The checks

After each question there is a `check_` call. It prints `OK` when your answer looks right,
or a hint when it does not. They are there to help you, not to grade you — nothing in this
project is scored or recorded. If a check disagrees with you and you think you are right,
you might be; read the hint and decide.

Written questions have checks too, but those only confirm you wrote something. Nobody is
marking your sentences.
"""

SOURCE_README = """\
# Regenerating the module

Everything in this bundle is generated. If the module needs editing, edit the source and
rebuild — do not edit the `.ipynb` files directly, because the next rebuild overwrites them.

| | |
|---|---|
| `build_dataset.py` | rebuilds `data/grains.csv` and the mystery files from the upstream FASTA downloads. Fixed random seed, so it reproduces exactly. |
| `build_notebooks.py` | the single source for all four notebooks. Both tracks and both versions (student / solutions) are generated from the same question text, so they cannot drift apart. |
| `rebuild.sh` | runs `build_notebooks.py`, then executes both solutions notebooks to prove they still work. Run this after any edit. |
| `package.py` | rebuilds the delivery zips. |

To change the wording of a question, find it in `build_notebooks.py` — the prose lives in
`md(r\"\"\"...\"\"\")` blocks and is plain Markdown — then run `./rebuild.sh`.

Requires the same packages as the notebooks, plus `jupyter` for the rebuild script.
"""


if __name__ == "__main__":
    raise SystemExit(main())
