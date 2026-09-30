# Data 8 Project — Spotting Gluten Grains from Protein Sequence

A self-contained data-science project built on a slice of a computational-biology
research project about celiac disease. Students never need biology background: a
protein is a string of letters, and the whole project is *counting letters, comparing
groups, and building a classifier*.

## The question

> Wheat, barley, and rye trigger celiac disease; rice, maize, oats, quinoa, millet, and
> buckwheat do not. The seed-storage proteins of the harmful grains are unusually rich in
> the amino acids **proline (P)** and **glutamine (Q)**. Can we use that — and a sharper
> sequence feature — to **predict** whether a protein comes from a gluten grain, and how
> well?

## Files

| File | What it is |
|---|---|
| `data/grains.csv` | 251 labeled seed-storage proteins: `seq_id, grain, gluten, protein_name, accession, sequence` |
| `data/grains_mystery.csv` | 18 **unlabeled** proteins students classify at the end |
| `data/grains_mystery_key.csv` | answers for the mystery set — **instructor only** |
| `build_dataset.py` | regenerates the three CSVs from the upstream FASTA files (reproducibility) |
| `INSTRUCTOR.md` | learning objectives, the five parts question-by-question, expected answers, review guidance, timing, common mistakes |
| `datascience_vs_pandas.md` | side-by-side reference for every table operation the project uses |
| `project_datascience.ipynb` | student notebook, Berkeley `datascience` library |
| `project_pandas.ipynb` | student notebook, `pandas` |
| `solutions_datascience.ipynb` / `solutions_pandas.ipynb` | full solutions |

## The five parts

0. **Orientation** — load the data, understand what a row is.
1. **Sequence features** — write functions that turn a sequence string into numbers
   (proline %, glutamine %, length).
2. **Exploration** — histograms and a scatter; which non-gluten grains look suspicious?
3. **Is the difference real?** — a permutation test on the gap in mean P+Q content.
4. **A sharper feature** — count a repeat motif; watch the borderline grains separate.
5. **Prediction** — a k-nearest-neighbors classifier: train/test split, accuracy vs. a
   baseline, feature-set comparison, error analysis, and the mystery set.

The emphasis is on **prediction** (Part 5). Judged on a held-out
"mystery set", because a random test split from the same grains is near-saturated.

## How to run it

```bash
python build_dataset.py        # only needed to rebuild data/ ; the CSVs are already here
jupyter notebook               # then open a project_*.ipynb
```

Needs `numpy`, `matplotlib`, and either `datascience` or `pandas` (all in the standard
Berkeley / Anaconda environment). No internet required.
