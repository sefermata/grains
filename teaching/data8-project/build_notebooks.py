"""Generate the project notebooks from a single source of truth.

Produces four .ipynb files, kept in sync:

    project_datascience.ipynb    student, Berkeley `datascience` library
    project_pandas.ipynb         student, pandas
    solutions_datascience.ipynb  full solutions, datascience
    solutions_pandas.ipynb       full solutions, pandas

Each content function returns a list of cells for one (track, version) pair:
    track   in {"ds", "pd"}
    version in {"project", "solutions"}      # "project" == scaffolded blanks

Authoring rules (so the nested-string situation stays sane):
  * every snippet below is a raw triple-double-quoted string
  * inside a snippet, use 'single quotes' for all Python string literals and f-strings
  * never put a triple-double-quote or an f"double" string inside a snippet

Run:  python build_notebooks.py
Currently emits Part 0 + Part 1 only (preview).
"""

from __future__ import annotations

import json
import textwrap
from pathlib import Path

HERE = Path(__file__).resolve().parent

NB_METADATA = {
    "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
    "language_info": {"name": "python", "version": "3"},
}


def _dedent(text: str) -> str:
    return textwrap.dedent(text).strip("\n") + "\n"


def md(text: str) -> dict:
    return {"cell_type": "markdown", "metadata": {}, "source": _dedent(text).splitlines(keepends=True)}


def code(text: str) -> dict:
    return {
        "cell_type": "code",
        "metadata": {},
        "execution_count": None,
        "outputs": [],
        "source": _dedent(text).splitlines(keepends=True),
    }


def solvable(version: str, *, scaffold: str, solution: str) -> dict:
    return code(solution if version == "solutions" else scaffold)


# --------------------------------------------------------------------------- #
#  front matter
# --------------------------------------------------------------------------- #

TITLE = r"""
# Project — Spotting Gluten Grains from Protein Sequence

Wheat, barley, and rye trigger celiac disease. Rice, maize, oats, quinoa, millet, and
buckwheat do not. In this project you'll test whether you can tell the two apart **from
protein sequence alone**, using only the tools from this course.
"""

OVERVIEW = r"""
## What you'll do

The project has six parts. Each one ends with something concrete you can check.

| Part | You will... | New skill |
|---|---|---|
| **0. Orientation** | Load the data and learn what one row means | reading a data table |
| **1. Sequence to numbers** | Write functions that turn a protein string into numbers (how much of it is proline, how much is glutamine) | functions, `apply` |
| **2. Exploration** | Compare gluten vs non-gluten grains with a histogram and a scatter plot; find which "safe" grains look suspicious | grouped visualization |
| **3. Is the difference real?** | Run a permutation test on the gap between the two groups | hypothesis testing |
| **5. A sharper feature** | Count a *repeating pattern* instead of raw amounts, and watch the borderline grains separate | more string logic |
| **6. Prediction** | Build a k-nearest-neighbors classifier, test it honestly on held-out data, and see where it fails | classification |

This preview contains **Part 0** and **Part 1**.
"""

WHY = r"""
## Why this project

**The real question.** Some grains harm people with celiac disease because of a protein
called **gluten**. Gluten is hard to digest, and fragments of it set off an immune reaction
in the small intestine. Researchers study exactly how those fragments fit into an immune
receptor using 3-D structure prediction. This project asks a simpler version of the same
question that fits in an intro course: *can the amino-acid makeup of a grain's proteins,
by itself, separate the harmful grains from the harmless ones?* And there's a built-in
answer key — we already know which grains are which — so you can tell whether your method
works.

**This is bioinformatics.** Bioinformatics (also called computational biology) means
answering biological questions with code and statistics applied to large datasets —
DNA, protein sequences, clinical measurements. The core skills are the ones this course
teaches: load and clean data, visualize it, test whether a pattern is real, and build a
model that predicts. The extra ingredient is a little biology vocabulary, which this
notebook gives you as you go.

**As a career.** Roles that combine data skills with biology have grown alongside the size
of biological datasets — in biotechnology and pharmaceutical companies, hospital and
public-health laboratories, agricultural research, and academic labs. You do not need a
biology degree to enter the field; a data-science foundation plus domain curiosity is a
common path. This project is a small proof of that: you'll reach a real result with nothing
beyond this course.
"""

BIOLOGY_MINIMUM = r"""
## The biology you need (all of it)

* A **protein** is a molecule built from a chain of small parts called **amino acids**.
  There are 20 kinds. Each is written as one capital letter:
  `A C D E F G H I K L M N P Q R S T V W Y`.
* So a protein can be written as a **string of letters**. For example, one real protein in
  our data begins:

  ```
  MAAKIFALLALHALLVSGTTAAIIPQCSLAPNAIIPQFIPPVTALGNEHLAVQAIRLQQ...
  ```

  That is not a misspelled word. It is a chain of 267 amino acids, one letter each, read
  left to right.
* Two of the letters matter here: **`P`** (the amino acid *proline*) and **`Q`**
  (*glutamine*). The storage protein in wheat, barley, and rye is unusually rich in `P`
  and `Q`, and that is chemically why it triggers celiac disease.
* In our data, each row is one protein. The `sequence` column holds its string. The
  `gluten` column says whether that protein's grain is one of the three that cause celiac
  disease. **Predicting `gluten` from `sequence` is the whole project.**
"""

TRACK_NOTE_DS = r"""
> **Track: `datascience`** — the Berkeley teaching library (`Table`). A parallel `pandas`
> version of this notebook exists; `datascience_vs_pandas.md` maps one to the other.
"""

TRACK_NOTE_PD = r"""
> **Track: `pandas`** — the library used outside the classroom. A parallel `datascience`
> version of this notebook exists; `datascience_vs_pandas.md` maps one to the other. The
> **pandas note** callouts flag places where pandas behaves differently from `datascience`.
"""

CHECKS_MD = r"""
## Setup — run these first

The next cell defines small `check_*` helpers used after every question. They print `OK`
when an answer looks right and a short hint otherwise. They are guidance, not a grade, and
you don't need to read them.
"""

SETUP_DS = r"""
from datascience import Table, are
import numpy as np

np.random.seed(8)   # keep every run reproducible (Data 8, naturally)
"""

SETUP_PD = r"""
import pandas as pd
import numpy as np

np.random.seed(8)   # keep every run reproducible (Data 8, naturally)
"""

# ---- checks that are pure Python (same in both tracks) -------------------- #
CHECKS_BASE = r"""
def _report(name, checks):
    for ok, hint in checks:
        if not ok:
            print('x  ' + name + ': ' + hint)
            return
    print('OK ' + name)


def check_written(name, text):
    bad = (not isinstance(text, str)) or len(text.strip()) < 15 or 'your answer' in text.lower()
    if bad:
        print('x  ' + name + ': write a sentence or two in the answer variable, replacing the placeholder')
    else:
        print('OK ' + name + '  (recorded - a person will read this one)')


def check_q0_1(n_rows, n_cols):
    _report('Q0.1', [
        (n_rows == 251, 'expected 251 rows, got ' + str(n_rows)),
        (n_cols == 6, 'expected 6 columns, got ' + str(n_cols)),
    ])


def check_q0_3(baseline_accuracy):
    ok = isinstance(baseline_accuracy, (int, float)) and 0.52 < baseline_accuracy < 0.54
    _report('Q0.3', [
        (ok, 'expected a fraction near 0.53 (133 of 251 rows), got ' + str(baseline_accuracy)),
    ])


def check_q1_1(fn):
    _report('Q1.1', [
        (fn('PQPQP', 'P') == 3, 'count_aa("PQPQP","P") should be 3, got ' + str(fn('PQPQP', 'P'))),
        (fn('PQPQP', 'Q') == 2, 'count_aa("PQPQP","Q") should be 2, got ' + str(fn('PQPQP', 'Q'))),
        (fn('AAAA', 'P') == 0, 'there is no P in "AAAA", so this should be 0'),
    ])


def check_q1_2(frac, pf, qf):
    try:
        got = [frac('PQPQ', 'P'), pf('PPPQ'), qf('PPPQ')]
    except Exception as e:
        print('x  Q1.2: a function raised an error - ' + str(e)); return
    if not all(isinstance(v, (int, float)) for v in got):
        print('x  Q1.2: each function should return a number (did you fill in the body?)'); return
    _report('Q1.2', [
        (abs(got[0] - 0.5) < 1e-9, 'aa_fraction("PQPQ","P") should be 0.5'),
        (abs(got[1] - 0.75) < 1e-9, 'proline_fraction("PPPQ") should be 0.75'),
        (abs(got[2] - 0.25) < 1e-9, 'glutamine_fraction("PPPQ") should be 0.25'),
    ])


def check_q1_4(longest_seq, shortest_seq, long_p_count, long_p_frac, short_p_count, short_p_frac):
    vals = [long_p_count, long_p_frac, short_p_count, short_p_frac]
    if not all(isinstance(v, (int, float)) for v in vals):
        print('x  Q1.4: set all four variables to numbers (use your Q1.1 / Q1.2 functions)'); return
    _report('Q1.4', [
        (long_p_count == longest_seq.count('P'), 'long_p_count should be the number of P in the longest sequence'),
        (short_p_count == shortest_seq.count('P'), 'short_p_count should be the number of P in the shortest sequence'),
        (abs(long_p_frac - longest_seq.count('P') / len(longest_seq)) < 1e-9, 'long_p_frac should be count / length'),
        (abs(short_p_frac - shortest_seq.count('P') / len(shortest_seq)) < 1e-9, 'short_p_frac should be count / length'),
    ])
"""

# ---- checks that touch the table / dataframe ---------------------------- #
CHECKS_DS_EXTRA = r"""
def check_q0_2(per_grain):
    try:
        n = per_grain.num_rows
        top = per_grain.sort('count', descending=True).column(0)[0]
    except Exception:
        print('x  Q0.2: per_grain should be a Table with one row per grain and a count column'); return
    _report('Q0.2', [
        (n == 10, 'expected 10 grains, got ' + str(n)),
        (top == 'wheat', 'the grain with the most proteins should be wheat, got ' + str(top)),
    ])


def check_q1_3(t):
    try:
        labels = set(t.labels)
    except Exception:
        print('x  Q1.3: run the column-adding code first'); return
    if not {'p_frac', 'q_frac', 'pq_frac'} <= labels:
        print('x  Q1.3: missing one of p_frac / q_frac / pq_frac'); return
    p, q, pq = t.column('p_frac'), t.column('q_frac'), t.column('pq_frac')
    gy = t.where('gluten', 'yes').column('pq_frac').mean()
    gn = t.where('gluten', 'no').column('pq_frac').mean()
    _report('Q1.3', [
        (p.min() >= 0 and p.max() <= 1, 'p_frac must be between 0 and 1'),
        (np.allclose(pq, p + q), 'pq_frac should equal p_frac + q_frac'),
        (0.40 < gy < 0.55, 'mean pq_frac for gluten grains looks off (' + format(gy, '.3f') + ', expect ~0.47)'),
        (gy > gn + 0.2, 'gluten grains should have a much higher mean pq_frac than the rest'),
    ])
"""

CHECKS_PD_EXTRA = r"""
def check_q0_2(per_grain):
    try:
        n = len(per_grain)
        top = per_grain.sort_values(ascending=False).index[0]
    except Exception:
        print('x  Q0.2: per_grain should be a Series of counts indexed by grain'); return
    _report('Q0.2', [
        (n == 10, 'expected 10 grains, got ' + str(n)),
        (top == 'wheat', 'the grain with the most proteins should be wheat, got ' + str(top)),
    ])


def check_q1_3(df):
    try:
        columns = set(df.columns)
    except Exception:
        print('x  Q1.3: run the column-adding code first'); return
    if not {'p_frac', 'q_frac', 'pq_frac'} <= columns:
        print('x  Q1.3: missing one of p_frac / q_frac / pq_frac'); return
    p, q, pq = df['p_frac'].to_numpy(), df['q_frac'].to_numpy(), df['pq_frac'].to_numpy()
    gy = df.loc[df['gluten'] == 'yes', 'pq_frac'].mean()
    gn = df.loc[df['gluten'] == 'no', 'pq_frac'].mean()
    _report('Q1.3', [
        (p.min() >= 0 and p.max() <= 1, 'p_frac must be between 0 and 1'),
        (np.allclose(pq, p + q), 'pq_frac should equal p_frac + q_frac'),
        (0.40 < gy < 0.55, 'mean pq_frac for gluten grains looks off (' + format(gy, '.3f') + ', expect ~0.47)'),
        (gy > gn + 0.2, 'gluten grains should have a much higher mean pq_frac than the rest'),
    ])
"""


# --------------------------------------------------------------------------- #
#  Part 0
# --------------------------------------------------------------------------- #

def part0(track: str, version: str) -> list[dict]:
    cells = [md("## Part 0 — Orientation  ·  5 points")]

    if track == "ds":
        cells.append(md("Load the data into a `Table` and look at the first few rows."))
        cells.append(code("grains = Table.read_table('data/grains.csv')\ngrains.show(8)"))
    else:
        cells.append(md("Load the data into a `DataFrame` and look at the first few rows."))
        cells.append(code("grains = pd.read_csv('data/grains.csv')\ngrains.head(8)"))

    cells.append(md(r"""
        ### What each column means

        | column | what it is | do we compute from it? |
        |---|---|---|
        | `seq_id` | our label for the row (`S001`, `S002`, ...) | No — it just lets us point at one specific protein. |
        | `grain` | the grain this protein comes from (`wheat`, `rice`, ...) | Used for grouping, and for error analysis in Part 6. |
        | `gluten` | `yes` if the grain is wheat, barley, or rye; otherwise `no` | **Yes — this is the label we try to predict.** |
        | `protein_name` | the name the source database gives this protein | No calculations from it. In Part 6 it helps explain *which* proteins the model gets wrong. Ignore it until then. |
        | `accession` | this protein's ID in **UniProt**, a public protein database | No — it's *provenance* only (a record of where the data came from, so you could trace it back). Paste it into uniprot.org to see the original record. |
        | `sequence` | the protein itself: a string of capital letters, one per amino acid | **Yes — every number we build comes from this column.** |

        You'll see the same `protein_name` (and the same `grain`) on many rows. That's
        expected: a grain has many storage-protein genes and many database entries, and
        each row is one protein record.
    """))

    cells.append(md(r"""
        ### Question 0.1 · 2 pt

        Set `n_rows` and `n_cols` to the number of rows and columns in `grains`.
    """))
    if track == "ds":
        cells.append(solvable(
            version,
            scaffold="n_rows = ...\nn_cols = ...\n\ncheck_q0_1(n_rows, n_cols)",
            solution="n_rows = grains.num_rows\nn_cols = grains.num_columns\n\ncheck_q0_1(n_rows, n_cols)",
        ))
    else:
        cells.append(solvable(
            version,
            scaffold="n_rows = ...\nn_cols = ...\n\ncheck_q0_1(n_rows, n_cols)",
            solution="n_rows = grains.shape[0]\nn_cols = grains.shape[1]\n\ncheck_q0_1(n_rows, n_cols)",
        ))

    cells.append(md(r"""
        ### Question 0.2 · 2 pt

        Each row is one protein, so counting rows per grain tells you how many proteins we
        have for each. Build a 2-column table, `per_grain`: one column is `grain`, the other
        is the count of proteins for that grain, sorted so the largest count comes first.
    """))
    if track == "ds":
        cells.append(solvable(
            version,
            scaffold="per_grain = ...\nper_grain",
            solution="per_grain = grains.group('grain').sort('count', descending=True)\nper_grain",
        ))
    else:
        cells.append(solvable(
            version,
            scaffold="per_grain = ...\nper_grain",
            solution="per_grain = grains['grain'].value_counts()\nper_grain",
        ))
        cells.append(md(r"""
            > **pandas note.** `value_counts()` returns a `Series` indexed by grain, already
            > sorted largest first. The `datascience` version returns a two-column `Table`
            > (`grain`, `count`). Same information, different shape.
        """))
    cells.append(code("check_q0_2(per_grain)"))

    cells.append(md(r"""
        ### Question 0.3 · 1 pt

        Before building any model, find the score to beat. Suppose you ignored the
        `sequence` entirely and, for every row, just guessed whichever `gluten` value shows
        up more often. How often would that be right?

        Count the `'no'` rows and the `'yes'` rows, then set `baseline_accuracy` to the
        larger count divided by the total number of rows.
    """))
    if track == "ds":
        cells.append(solvable(
            version,
            scaffold=r"""
                n_no  = ...   # number of rows where gluten == 'no'
                n_yes = ...   # number of rows where gluten == 'yes'
                baseline_accuracy = ...   # bigger of the two counts, divided by the total number of rows

                print('always guess the majority label ->', baseline_accuracy)
                check_q0_3(baseline_accuracy)
            """,
            solution=r"""
                n_no  = grains.where('gluten', 'no').num_rows
                n_yes = grains.where('gluten', 'yes').num_rows
                baseline_accuracy = max(n_no, n_yes) / grains.num_rows

                print('always guess the majority label ->', baseline_accuracy)
                check_q0_3(baseline_accuracy)
            """,
        ))
    else:
        cells.append(solvable(
            version,
            scaffold=r"""
                n_no  = ...   # number of rows where gluten == 'no'
                n_yes = ...   # number of rows where gluten == 'yes'
                baseline_accuracy = ...   # bigger of the two counts, divided by the total number of rows

                print('always guess the majority label ->', baseline_accuracy)
                check_q0_3(baseline_accuracy)
            """,
            solution=r"""
                n_no  = (grains['gluten'] == 'no').sum()
                n_yes = (grains['gluten'] == 'yes').sum()
                baseline_accuracy = max(n_no, n_yes) / len(grains)

                print('always guess the majority label ->', baseline_accuracy)
                check_q0_3(baseline_accuracy)
            """,
        ))
    cells.append(md(r"""
        Now say why this number matters, in one or two sentences. (What would it mean if a
        model you build only reached this accuracy?) Put your answer in `answer_q0_3`.
    """))
    cells.append(solvable(
        version,
        scaffold="answer_q0_3 = 'your answer here'\n\ncheck_written('Q0.3 written', answer_q0_3)",
        solution=(
            "answer_q0_3 = (\n"
            "    'It is the accuracy of a model that ignores the sequence completely, so any real '\n"
            "    'model has to beat it to be worth anything. Reaching only this number would mean '\n"
            "    'the model learned the class sizes, not the biology.'\n"
            ")\n\ncheck_written('Q0.3 written', answer_q0_3)"
        ),
    ))
    return cells


# --------------------------------------------------------------------------- #
#  Part 1
# --------------------------------------------------------------------------- #

def part1(track: str, version: str) -> list[dict]:
    cells = [md(r"""
        ## Part 1 — From a sequence to numbers  ·  12 points

        A model can't work with the raw string `'MAAKIF...IF'`. It needs **features**:
        numbers that summarize a sequence. In this part you build three of them —

        * `p_frac`: the fraction of the protein that is proline (`P`)
        * `q_frac`: the fraction that is glutamine (`Q`)
        * `pq_frac`: the two added together

        — and attach them to the table as new columns.

        You'll start with a general helper, `count_aa(seq, aa)`. The name `aa` is short for
        **amino acid**: it's just a parameter, so the one function works for any letter. In
        Q1.2 you'll call it with `'P'` and `'Q'` specifically. We don't care about amino
        acids in general — `aa` is only there so we don't have to write the counting logic
        twice.
    """)]

    # Q1.1
    cells.append(md(r"""
        ### Question 1.1 · 3 pt

        Write `count_aa(seq, aa)` that returns how many times the single character `aa`
        occurs in the string `seq`. A `for` loop over the letters of `seq` is a fine way to
        do it.
    """))
    cells.append(solvable(
        version,
        scaffold=r"""
            def count_aa(seq, aa):
                # number of times the character aa appears in the string seq
                ...

            check_q1_1(count_aa)
        """,
        solution=r"""
            def count_aa(seq, aa):
                # number of times the character aa appears in the string seq
                aa_count = 0
                for letter in seq:
                    if letter == aa:
                        aa_count = aa_count + 1
                return aa_count

            # (a one-line version, once the loop feels obvious: sum(letter == aa for letter in seq))
            check_q1_1(count_aa)
        """,
    ))

    # Q1.2
    cells.append(md(r"""
        ### Question 1.2 · 3 pt

        Write `aa_fraction(seq, aa)`: the fraction of `seq` that is the character `aa` (a
        number between 0 and 1). Reuse `count_aa` rather than re-writing the counting logic.

        Then write `proline_fraction(seq)` and `glutamine_fraction(seq)` by calling
        `aa_fraction` with the right letter for each amino acid.

        > **Note:** string comparisons in Python are case-sensitive — `'P'` and `'p'` are
        > different characters, so double-check which case the `sequence` column actually uses.
    """))
    cells.append(solvable(
        version,
        scaffold=r"""
            def aa_fraction(seq, aa):
                ...

            def proline_fraction(seq):
                ...

            def glutamine_fraction(seq):
                ...

            check_q1_2(aa_fraction, proline_fraction, glutamine_fraction)
        """,
        solution=r"""
            def aa_fraction(seq, aa):
                return count_aa(seq, aa) / len(seq)

            def proline_fraction(seq):
                return aa_fraction(seq, 'P')

            def glutamine_fraction(seq):
                return aa_fraction(seq, 'Q')

            check_q1_2(aa_fraction, proline_fraction, glutamine_fraction)
        """,
    ))

    # re-show the table before adding columns
    cells.append(md("Here is the table again, so the column names are fresh before you add to it:"))
    cells.append(code("grains.show(3)" if track == "ds" else "grains.head(3)"))

    # Q1.3
    if track == "ds":
        cells.append(md(r"""
            ### Question 1.3 · 4 pt

            Add three columns to `grains`: `p_frac`, `q_frac`, and `pq_frac`
            (= `p_frac + q_frac`), computed from the `sequence` column with the functions
            you just wrote.

            **Hints:** `with_columns` returns a *new* table rather than changing `grains` in
            place — make sure the result actually gets saved back to `grains`. A `Table` also
            has an `apply` method that runs a function down every value in one column; look
            it up in the `datascience` docs if you haven't used it yet.
        """))
        cells.append(solvable(
            version,
            scaffold=r"""
                grains = grains.with_columns(
                    'p_frac', ...,
                    'q_frac', ...,
                )
                grains = grains.with_columns('pq_frac', ...)

                grains.select('seq_id', 'grain', 'gluten', 'p_frac', 'q_frac', 'pq_frac').show(5)
                check_q1_3(grains)
            """,
            solution=r"""
                grains = grains.with_columns(
                    'p_frac', grains.apply(proline_fraction, 'sequence'),
                    'q_frac', grains.apply(glutamine_fraction, 'sequence'),
                )
                grains = grains.with_columns('pq_frac', grains.column('p_frac') + grains.column('q_frac'))

                grains.select('seq_id', 'grain', 'gluten', 'p_frac', 'q_frac', 'pq_frac').show(5)
                check_q1_3(grains)
            """,
        ))
    else:
        cells.append(md(r"""
            ### Question 1.3 · 4 pt

            Add three columns to `grains`: `p_frac`, `q_frac`, and `pq_frac`
            (= `p_frac + q_frac`), computed from the `sequence` column with the functions
            you just wrote.

            **Hints:** a pandas `Series` (a single column) has an `apply` method that runs a
            function on every value in it. Assigning to `grains['a_new_column']` adds that
            column in place — no reassigning `grains` needed here, unlike the `datascience`
            track.
        """))
        cells.append(solvable(
            version,
            scaffold=r"""
                grains['p_frac'] = ...
                grains['q_frac'] = ...
                grains['pq_frac'] = ...

                grains[['seq_id', 'grain', 'gluten', 'p_frac', 'q_frac', 'pq_frac']].head(5)
                check_q1_3(grains)
            """,
            solution=r"""
                grains['p_frac'] = grains['sequence'].apply(proline_fraction)
                grains['q_frac'] = grains['sequence'].apply(glutamine_fraction)
                grains['pq_frac'] = grains['p_frac'] + grains['q_frac']

                grains[['seq_id', 'grain', 'gluten', 'p_frac', 'q_frac', 'pq_frac']].head(5)
                check_q1_3(grains)
            """,
        ))

    # Q1.4 — concrete, then a one-sentence written follow-up
    seq_access = "grains.column('sequence')" if track == "ds" else "grains['sequence']"
    cells.append(md(r"""
        ### Question 1.4 · 2 pt

        Here's why we used a *fraction* and not a raw count. The two cells below pick out the
        longest and the shortest protein in the dataset. For each, set the raw `P` count and
        the proline fraction, using your functions from Q1.1 and Q1.2.
    """))
    cells.append(code(f"longest_seq = max({seq_access}, key=len)\nshortest_seq = min({seq_access}, key=len)\nprint('lengths:', len(longest_seq), 'vs', len(shortest_seq))"))
    cells.append(solvable(
        version,
        scaffold=r"""
            long_p_count  = ...
            long_p_frac   = ...
            short_p_count = ...
            short_p_frac  = ...

            print('longest :', long_p_count, 'P out of', len(longest_seq), '->', long_p_frac)
            print('shortest:', short_p_count, 'P out of', len(shortest_seq), '->', short_p_frac)
            check_q1_4(longest_seq, shortest_seq, long_p_count, long_p_frac, short_p_count, short_p_frac)
        """,
        solution=r"""
            long_p_count  = count_aa(longest_seq, 'P')
            long_p_frac   = proline_fraction(longest_seq)
            short_p_count = count_aa(shortest_seq, 'P')
            short_p_frac  = proline_fraction(shortest_seq)

            print('longest :', long_p_count, 'P out of', len(longest_seq), '->', long_p_frac)
            print('shortest:', short_p_count, 'P out of', len(shortest_seq), '->', short_p_frac)
            check_q1_4(longest_seq, shortest_seq, long_p_count, long_p_frac, short_p_count, short_p_frac)
        """,
    ))
    cells.append(md(r"""
        You want to say which of these two proteins is *more proline-rich*. The raw counts
        and the fractions can point to different answers. Which number answers that question
        fairly, and why? One sentence in `answer_q1_4`.
    """))
    cells.append(solvable(
        version,
        scaffold="answer_q1_4 = 'your answer here'\n\ncheck_written('Q1.4 written', answer_q1_4)",
        solution=(
            "answer_q1_4 = (\n"
            "    'The fraction, because it divides out length: a raw count is large for any long '\n"
            "    'protein and small for any short one, so it mostly measures size, not how P-rich '\n"
            "    'the protein is.'\n"
            ")\n\ncheck_written('Q1.4 written', answer_q1_4)"
        ),
    ))

    cells.append(md(r"""
        ---
        **Part 1 done.** `grains` now carries `p_frac`, `q_frac`, and `pq_frac`. Part 2
        uses them to compare gluten vs non-gluten grains with a histogram and a scatter
        plot — and to spot which "safe" grains look suspiciously gluten-like.
    """))
    return cells


# --------------------------------------------------------------------------- #
#  assembly
# --------------------------------------------------------------------------- #

def build(track: str, version: str) -> dict:
    checks = CHECKS_BASE + "\n" + (CHECKS_DS_EXTRA if track == "ds" else CHECKS_PD_EXTRA)
    cells = [
        md(TITLE),
        md(OVERVIEW),
        md(WHY),
        md(BIOLOGY_MINIMUM),
        md(TRACK_NOTE_DS if track == "ds" else TRACK_NOTE_PD),
        md(CHECKS_MD),
        code(SETUP_DS if track == "ds" else SETUP_PD),
        code(checks),
        *part0(track, version),
        *part1(track, version),
    ]
    return {"cells": cells, "metadata": NB_METADATA, "nbformat": 4, "nbformat_minor": 4}


def main() -> None:
    targets = {
        ("ds", "project"): "project_datascience.ipynb",
        ("pd", "project"): "project_pandas.ipynb",
        ("ds", "solutions"): "solutions_datascience.ipynb",
        ("pd", "solutions"): "solutions_pandas.ipynb",
    }
    for (track, version), fname in targets.items():
        nb = build(track, version)
        (HERE / fname).write_text(json.dumps(nb, indent=1))
        print(f"wrote {fname}  ({len(nb['cells'])} cells)")


if __name__ == "__main__":
    main()
