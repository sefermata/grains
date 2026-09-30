# `datascience` vs `pandas` — the operations this project uses

The Berkeley `datascience` library was built for teaching: one data structure (`Table`),
verbose method names, few surprises. `pandas` is what practitioners use: more powerful, more
concepts (`Index`, `Series`, `dtype`, views vs copies), more ways to do the same thing.

This project is written in both. Same analysis, same numbers, same figures. Below is every
operation the notebooks use, side by side.

Setup:

```python
# datascience track
from datascience import Table
import numpy as np

# pandas track
import pandas as pd
import numpy as np
```

---

## Load a CSV

```python
t = Table.read_table("data/grains.csv")            # datascience
df = pd.read_csv("data/grains.csv")                # pandas
```

## Shape / peek

| task | datascience | pandas |
|---|---|---|
| first rows | `t.show(3)` | `df.head(3)` |
| row count | `t.num_rows` | `len(df)` / `df.shape[0]` |
| column count | `t.num_columns` | `df.shape[1]` |
| column names | `t.labels` | `df.columns` |
| one column as an array | `t.column("grain")` | `df["grain"].to_numpy()` |

## Filter rows

```python
# rows where gluten == "yes"
t.where("gluten", "yes")                           # datascience
df[df["gluten"] == "yes"]                          # pandas

# rows where pq_frac > 0.3
t.where("pq_frac", are.above(0.3))                 # datascience  (from datascience import are)
df[df["pq_frac"] > 0.3]                            # pandas
```

## Add / compute a column

```python
# datascience — returns a NEW table, must reassign
t = t.with_columns("p_frac", t.apply(proline_fraction, "sequence"))

# pandas — several ways; this one is closest in spirit
df["p_frac"] = df["sequence"].apply(proline_fraction)
# or, functional style that also returns a new frame:
df = df.assign(p_frac=df["sequence"].apply(proline_fraction))
```

`Table.apply(f, "col")` calls `f` on each value of one column and returns an array.
`Series.apply(f)` is the direct equivalent. For a function of several columns:

```python
t.apply(lambda seq, n: some_f(seq, n), "sequence", "n_residues")     # datascience
df.apply(lambda row: some_f(row["sequence"], row["n_residues"]), axis=1)   # pandas
```

## Group and aggregate

```python
# count rows per grain
t.group("grain")                                   # datascience -> columns: grain, count
df["grain"].value_counts()                         # pandas -> a Series
df.groupby("grain").size()                         # pandas -> a Series

# mean of pq_frac by gluten
t.select("gluten", "pq_frac").group("gluten", np.mean)      # datascience
df.groupby("gluten")["pq_frac"].mean()                       # pandas
```

## Sort

```python
t.sort("count", descending=True)                   # datascience
df.sort_values("count", ascending=False)           # pandas
```

## Select / drop columns

```python
t.select("grain", "sequence")                      # datascience
df[["grain", "sequence"]]                          # pandas

t.drop("accession")                                # datascience
df.drop(columns="accession")                       # pandas
```

## Sample / shuffle / split

```python
# shuffle all rows
shuf = t.sample(t.num_rows, with_replacement=False)          # datascience
shuf = df.sample(frac=1, random_state=8)                     # pandas

# 70/30 split of an already-shuffled table
n = int(0.7 * shuf.num_rows)
train, test = shuf.take(np.arange(n)), shuf.take(np.arange(n, shuf.num_rows))   # datascience
train, test = shuf.iloc[:n], shuf.iloc[n:]                                       # pandas

# bootstrap resample
t.sample(t.num_rows, with_replacement=True)                  # datascience
df.sample(n=len(df), replace=True, random_state=...)         # pandas
```

## Histogram / scatter

```python
# overlaid histogram of pq_frac by group
t.hist("pq_frac", group="gluten", bins=np.arange(0, 0.76, 0.02))      # datascience

# pandas: use matplotlib directly
import matplotlib.pyplot as plt
for label, sub in df.groupby("gluten"):
    plt.hist(sub["pq_frac"], bins=np.arange(0, 0.76, 0.02), alpha=0.5, label=label)
plt.legend(); plt.xlabel("pq_frac");

# scatter colored by group
t.scatter("p_frac", "q_frac", group="gluten")                        # datascience
for label, sub in df.groupby("gluten"):
    plt.scatter(sub["p_frac"], sub["q_frac"], s=10, alpha=0.6, label=label)
plt.legend(); plt.xlabel("p_frac"); plt.ylabel("q_frac");
```

## Row access inside a classifier

```python
# datascience: a row behaves like a named tuple
row = test.row(0)
row.item("p_frac")

# pandas: iterate rows as Series
row = test.iloc[0]
row["p_frac"]
# faster for many rows: work with numpy arrays directly
X_train = train[["p_frac", "q_frac"]].to_numpy()
```

## k-NN — hand-built vs library

The `datascience` track builds k-NN from the definition (distance function, majority vote),
matching how Data 8 teaches it. The `pandas` track does the same **and then** shows the
practitioner version:

```python
from sklearn.neighbors import KNeighborsClassifier

features = ["pq_motif", "p_frac"]
knn = KNeighborsClassifier(n_neighbors=5)
knn.fit(train[features].to_numpy(), train["gluten"].to_numpy())
pred = knn.predict(test[features].to_numpy())
accuracy = (pred == test["gluten"].to_numpy()).mean()
```

Same idea, ~5 lines, and it scales. The point of building it by hand first is to know
exactly what those 5 lines do.

---

## Gotchas when moving from `datascience` to `pandas`

| | |
|---|---|
| **Return vs mutate** | `Table` methods almost always return a new table. `pandas` has both — `df.sort_values(...)` returns a copy, `df.sort_values(..., inplace=True)` mutates. Prefer the copy-returning form while learning. |
| **The Index** | `pandas` attaches a row index; after filtering, `df.iloc[0]` (position) and `df.loc[0]` (label) can differ. `reset_index(drop=True)` after a filter avoids surprises. |
| **Series vs array** | `df["x"]` is a `Series` (has its own index); `df["x"].to_numpy()` is a plain array like `t.column("x")`. Mixing Series with mismatched indices aligns on the index, which can silently produce `NaN`. |
| **`SettingWithCopyWarning`** | assigning into a slice (`df[df.a>0]["b"] = ...`) may not stick. Use `.loc[mask, "b"] = ...` or `df = df.assign(...)`. |
| **Chained `groupby`** | `datascience` `group` always returns a table; `pandas` `groupby` returns a lazy object you must aggregate (`.mean()`, `.size()`, `.agg(...)`). |
