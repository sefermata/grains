# Instructor Guide — Spotting Gluten Grains from Protein Sequence

## Summary

Students build, from scratch, a pipeline that predicts whether a seed-storage protein comes
from a gluten grain (wheat / barley / rye) using only its amino-acid sequence. The project
walks from "turn a string into features" through a permutation test and a bootstrap interval
to a hand-built k-nearest-neighbors classifier with error analysis and a held-out mystery set.

It is a scoped-down version of a real research project (structural modeling of gluten
peptides against the HLA-DQ2.5 immune receptor). Students see the "just count the amino
acids" version of a question that researchers attack with 3D structure prediction.

**No biology background required.** Everything students need is stated in Part 0.

## Learning objectives

By the end, a student can:

1. Write pure-Python functions that compute summary features from a string, and apply them
   across a whole table (`Table.apply` / `df.apply`).
2. Compare two groups with overlaid histograms, a scatter, and group means, and describe the
   overlap in words.
3. Run a permutation test: state the null, shuffle labels, build the null distribution,
   compute an empirical p-value, and interpret it correctly (including what it does *not* say).
4. Build a bootstrap confidence interval and explain how it differs in purpose from the test.
5. Explain why a repeat-motif rate can be a more *specific* feature than raw composition,
   and demonstrate it on the borderline grains.
6. Implement k-NN classification from its definition: a distance function, a classifier that
   takes a majority vote of the `k` nearest training rows, a train/test split, and accuracy.
7. Compare classifiers against a majority-class baseline, compare feature sets on held-out
   accuracy, do error analysis by group, and apply the final model to new unlabeled data.
8. (pandas track) recognize the same operations in `datascience` and `pandas`, and use
   `scikit-learn`'s `KNeighborsClassifier` as the "in practice" version of the hand-built one.

## Format

Between fill-in-the-blank and from-scratch — modeled on the Data 8 projects (Movie
Classification, etc.), not the weekly labs. Each part has a short narrative, numbered
questions, some starter code, and written-answer cells. Roughly 60% code the student writes
with a function signature or one line of scaffolding given, 40% they structure themselves.

Lightweight self-checks are inline (`check_q1_2()` prints ✓ / a hint) so the notebook works
without `otter-grader`. To convert to otter: move each `check_*` body into a `tests/qN_M.py`
file and replace the calls with `grader.check("qN_M")`.

## Timing

| Use | Parts | Time |
|---|---|---|
| One long lab | 0–3 | 2–3 h |
| Checkpoint (due mid-project) | 0–4 | first week |
| Full project | 0–6 | 8–12 h over ~2 weeks |

## Dataset notes

- `grains.csv`: 251 rows, **118 gluten / 133 gluten-free** → majority-class baseline **53%**.
  Built with a fixed seed (`build_dataset.py`, `SEED = 8`); re-running reproduces it exactly.
- Grain-level mean P+Q fraction (what students will rediscover):

  | grain | n | P+Q frac | "PQ" per 100 aa |
  |---|--:|--:|--:|
  | rye | 3 | 0.54 | 11.4 |
  | barley | 35 | 0.47 | 7.9 |
  | wheat | 80 | 0.47 | 5.0 |
  | **oat** | 3 | **0.30** | **1.5** |
  | **maize** | 25 | **0.27** | 0.8 |
  | **sorghum** | 20 | **0.26** | 0.8 |
  | foxtail_millet | 20 | 0.14 | 0.2 |
  | quinoa | 20 | 0.13 | 0.3 |
  | buckwheat | 20 | 0.13 | 0.4 |
  | rice | 25 | 0.12 | 0.2 |

- **The teaching tension:** maize, sorghum, and oat have P/Q-rich storage proteins too
  (P+Q ≈ 0.27–0.30 vs ≈ 0.47 gluten, ≈ 0.13 other safe). On a random test split the classes
  still separate cleanly, so composition looks perfect — the cost only shows on the harder
  **mystery set**. The **"PQ" motif rate** (Part 5) pulls maize/sorghum down to the safe
  cluster (≈ 0.8 vs ≈ 6.0 for gluten) because their P and Q are not in the tandem
  `...PQ PQ PQ...` repeats characteristic of true gluten. Oat avenin stays mildly elevated —
  a genuine biological borderline, not a bug.
- `rye` has only 3 sequences (public databases under-annotate the rye genome). Enough that the
  gluten class is not purely wheat; small enough to mention as a real-world data limitation.
- Sequences are cleaned to the 20 standard amino acids, length 50–600, de-duplicated.
- Mystery set (18): 10 gluten (wheat LMW/HMW glutenin, gliadin, B/D/B3 hordein) + 8 gluten-free
  (rice glutelin, alpha-zein, and **two oat avenins, M07 & M18** — the intended "hmm" cases).

---

# The six parts, question by question

Point values in **[brackets]**. Total **105**, prediction-heavy by design: Parts 5–6 = 55.
(Part 4 is the trim point if you need to get back to 100 or shorten the project.)

## Part 0 — Orientation  **[5]**

Narrative gives students these facts: proteins are strings over a 20-letter amino-acid
alphabet; seeds store protein to feed the sprout; in wheat/barley/rye a large share of that
storage protein is gluten, which is proline- and glutamine-rich and hard to digest; that
same chemistry is what makes gluten reactive in celiac disease; the goal is to predict the
`gluten` label from `sequence`.

- **Q0.1 [2]** Load `data/grains.csv`. How many rows? How many columns?
  *Answer: 251 rows, 6 columns.*
- **Q0.2 [2]** Make a table of the number of proteins per `grain`, sorted descending.
  *Answer: wheat 80, barley 35, rice 25, maize 25, … rye 3, oat 3.*
- **Q0.3 [1]** *(written)* If a classifier ignored the sequence entirely and always guessed
  the more common label, what accuracy would it get? *Answer: 133/251 ≈ 53%. This is the
  baseline every later model must beat.*

## Part 1 — Sequence features  **[12]**

- **Q1.1 [3]** Write `count_aa(seq, aa)` → number of times character `aa` appears in `seq`.
- **Q1.2 [3]** Write `aa_fraction(seq, aa)` → `count_aa / len(seq)`. Use it to define
  `proline_fraction(seq)` and `glutamine_fraction(seq)`.
- **Q1.3 [4]** Add three columns to the table: `p_frac`, `q_frac`, `pq_frac`
  (= `p_frac + q_frac`), via `apply`.
- **Q1.4 [2]** *(written)* Why use a *fraction* rather than a raw count? *Answer: sequences
  vary in length (50–600 here); raw counts would mostly measure protein size, not composition.*

Common mistakes: using Python's `str.count` is fine here (non-overlapping is OK for single
characters); forgetting to divide by length; adding columns without reassigning the table in
the `datascience` track (`t = t.with_columns(...)`).

## Part 2 — Exploration  **[13]**

- **Q2.1 [4]** Overlaid histogram of `pq_frac`, one curve per `gluten` value. Reasonable bins.
- **Q2.2 [3]** Scatter `p_frac` (x) vs `q_frac` (y), colored by `gluten`.
- **Q2.3 [3]** Table of mean `p_frac`, `q_frac`, `pq_frac` grouped by `gluten`.
  *Answer: gluten yes ≈ 0.47; no ≈ 0.18.*
- **Q2.4 [3]** *(written)* Name the non-gluten grains whose points sit closest to the gluten
  cluster, and say why that matters for a classifier. *Answer: maize and sorghum (and oat).
  Their storage proteins (zein, kafirin, avenin) are also prolamins — P/Q-rich — so
  composition alone will produce false positives.*

## Part 3 — Is the difference real?  **[12]**

- **Q3.1 [2]** Compute the observed difference `mean(pq_frac | gluten) − mean(pq_frac | not)`.
- **Q3.2 [2]** *(written)* State the null hypothesis for a permutation test here.
  *Answer: the `gluten` label is unrelated to `pq_frac`; the observed gap is what you'd get
  from a random split of these same numbers into groups of 118 and 133.*
- **Q3.3 [5]** Write a simulation: shuffle the `gluten` column, recompute the difference in
  means, repeat ≥ 5000 times, collect the results.
- **Q3.4 [2]** Histogram of the simulated differences with the observed value marked; compute
  the empirical p-value. *Answer: no simulated difference comes near the observed ≈ 0.29;
  p ≈ 0 (report as "< 1/5000").*
- **Q3.5 [1]** *(written)* What does this p-value let you conclude, and what does it **not**
  tell you? *Answer: the composition gap is not plausibly due to chance labeling; it does
  **not** tell you the effect is large, nor that composition will classify well.*

Common mistakes: shuffling `pq_frac` instead of the label (equivalent here, but muddles the
logic); recomputing against the *shuffled* group sizes wrong; using `np.random.shuffle`
in place on a shared array without copying.

## Part 4 — How big is the difference?  **[8]**

Deliberately short — the course already drills bootstrap.

- **Q4.1 [4]** Bootstrap the difference in mean `pq_frac` (resample rows with replacement,
  ≥ 5000 times).
- **Q4.2 [2]** Report a 95% confidence interval from the 2.5th / 97.5th percentiles.
- **Q4.3 [2]** *(written)* One sentence: how is the question in Part 4 different from Part 3?
  *Answer: Part 3 asks "could this be nothing?"; Part 4 asks "given it's something, how large,
  with what uncertainty?"*

## Part 5 — A sharper feature  **[15]**

- **Q5.1 [5]** Write `motif_rate(seq, motif)` = number of **overlapping** occurrences of
  `motif` per 100 residues. Students must handle overlap themselves
  (`"PQPQ"` contains `"PQ"` twice; `str.count` gets non-overlapping only — for `"PQ"` that
  happens to match, so **test with `motif="QQ"` on `"QQQ"`**, which must give 2, not 1).
- **Q5.2 [3]** Add a `pq_motif` column = `motif_rate(seq, "PQ")`.
- **Q5.3 [3]** Group means of `pq_motif` by `gluten`, and by `grain`.
  *Answer: gluten ≈ 6.0 vs 0.5 per 100; maize and sorghum drop to ≈ 0.8 — now with the safe
  grains, not between.*
- **Q5.4 [4]** *(written)* Why is a repeat-dipeptide rate more *specific* to true gluten than
  the P+Q fraction? *Answer: zein/kafirin are P/Q-rich but their residues are not in tandem
  `PQ` repeats; the immunogenic gluten peptides come precisely from long `...PQ-PQ-PQ...`
  repeat tracts, so the motif rate captures arrangement, not just amount.*

## Part 6 — Prediction  **[35]**

> Reality check baked into this part: on a random test split from the **same grains**, a
> composition classifier is already near-perfect, so test accuracy barely separates good
> features from bad. The **mystery set** — deliberately made of atypical proteins (short
> HMW/D glutenins, oat avenins) — is where feature choice and feature scaling actually show.
> Students should come away trusting the mystery-set number over the test-set number.

### 6a. Setup

- **Q6.1 [3]** Shuffle the table (seed given) and split 70/30 into `train` / `test`.
- **Q6.2 [3]** *(written)* Why must the split happen before you look at test accuracy, and why
  shuffle first? *Answer: to estimate performance on unseen data; the CSV is grain-blocked, so
  an unshuffled split would put whole grains entirely on one side.*

### 6b. k-NN from scratch

- **Q6.3 [4]** `distance(row1, row2, features)` — Euclidean distance over a list of feature
  columns.
- **Q6.4 [5]** `classify(test_row, train, features, k)` — distance to every training row, take
  the `k` closest, return the majority `gluten` label.
- **Q6.5 [3]** `accuracy(train, evalset, features, k)`. Report it on `test` for
  `features = ["p_frac", "q_frac"]`, `k = 5`.
  *Answer: ≈ 0.98–1.00 depending on split — far above the 53% baseline, and suspiciously high.*
- **Q6.6 [2]** *(written)* Your test accuracy is ~99%. Give one reason not to trust that as the
  accuracy you'd get on protein from a grain not in this dataset. *Answer: train and test are
  random draws from the same 10 grains with the same protein families; the model can lean on
  grain-specific regularities that won't transfer.*

### 6c. Features and scaling — judged on the mystery set

- **Q6.7 [3]** Load `grains_mystery.csv`; compute the same feature columns on it. Add an
  `n_residues` column to `train`, `test`, and `mystery`.
- **Q6.8 [4]** For `k = 5`, fill a small table of **test** and **mystery** accuracy for:
  `["n_residues"]`, `["p_frac", "q_frac"]`, `["pq_motif", "p_frac"]` (raw).
  *Answer (≈): length 0.87 / 0.67; composition 0.99 / 0.83; motif+comp raw 0.92 / 0.83. Test
  accuracy hides the differences; the mystery column shows length fails to generalise.*
- **Q6.9 [4]** standardize every feature (subtract the **train** mean, divide by the **train**
  SD — the same transform applied to test and mystery). Re-evaluate `["pq_motif", "p_frac"]`.
  *Answer: mystery ≈ 0.83 → ≈ 0.89. Raw, `pq_motif` (range ~0–12) swamps `p_frac` (range
  ~0–0.2) in the distance; standardizing gives them equal say.*
- **Q6.10 [3]** *(written)* From these numbers, which feature set + scaling would you ship, and
  on which accuracy figure did you base that? *Answer: standardized `["pq_motif", "p_frac"]`,
  chosen on the mystery-set accuracy, since the test figure is saturated.*

### 6d. Error analysis and the reveal

- **Q6.11 [2]** With your chosen model, list the misclassified **mystery** proteins with their
  true `grain` (instructor reveals `grains_mystery_key.csv`).
  *Answer: usually 2 misses out of 18 — an oat avenin (M07/M18) called "gluten", and/or a
  short HMW/D-type glutenin or hordein (M01/M10) called "gluten-free".*
- **Q6.12 [2]** *(written)* Pick one miss and explain what about its sequence caused it.
  *Answer: oat avenin really does carry `PQ`-rich repeats (avenin is a prolamin) — the model
  is "wrong" about a case biology also treats as a gray area. HMW/D glutenins are gluten but
  glycine-rich with fewer tight `PQ` repeats, so they look composition-wise like a safe grain.*

### 6e. Open-ended close

- **Q6.13 [2]** *(written, ~150 words)* Propose one additional feature to raise mystery-set
  accuracy, name the grains it should help and why, and say how you'd check it genuinely helps
  rather than just sounding plausible. Strong answers: a concrete feature (rate of `"QPQ"` or
  a length-≥6 `PQ`-repeat run; fraction of residues inside any run of ≥4 alternating P/Q;
  presence of a known epitope core such as `"PFPQPQLPY"`), a mechanism linking it to gluten,
  and a real evaluation (re-check on the mystery set or cross-validate on `train` — not
  re-scoring the data used to pick it).

---

## Rubric (per written question)

| | Full | Partial | None |
|---|---|---|---|
| Written | Correct claim **and** correct reason; addresses what was asked | Right idea, vague or missing the "why" | Wrong, or restates the question |
| Code | Correct output **and** matches the required signature/approach | Right approach, an off-by-one / wrong axis / needs a nudge | Doesn't run or wrong method |

Cap Part 6 at full marks if 6c–6d work even when 6b is imperfect but reasonable — the point is
the prediction workflow, not a flawless distance function.

## Common student mistakes (collected)

- Adding a column in the `datascience` track without reassigning (`t.with_columns` returns a
  new table).
- `str.count` for overlapping motifs (Q5.1) — the deliberate trap; test on `"QQQ"`.
- Permutation test: shuffling a feature not the label; not enough repetitions; comparing to a
  one-sided vs two-sided threshold without saying which.
- Train/test leakage: computing `pq_motif` or picking `k` after peeking at `test`.
- Reporting accuracy without the 53% baseline next to it.
- Trusting the ~99% test accuracy — the point of Q6.6 / the mystery set is that it doesn't
  generalise.
- Q6.9 scaling: computing the mean/SD from the *whole* dataset or from `test`/`mystery`
  instead of from `train` only (a leak); forgetting to apply the same transform to
  `test`/`mystery`; standardizing the label column by accident.
- Euclidean distance left on raw features of different scales (`pq_motif` ~0–12 vs `p_frac`
  ~0–0.2) so the classifier effectively only sees `pq_motif`.
