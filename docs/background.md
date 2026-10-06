# Background: what this project is and why

Orientation written to be read cold. If you've been away from the repo for three weeks, start
here. Most of this becomes the introduction section of the final report.

Terms are defined where they first appear, and collected again in the [glossary](#glossary) at
the end.

---

## 1. What is physically happening in celiac disease

Someone eats bread. Their immune system attacks their own small intestine.

The lining of the small intestine is covered in finger-like projections called **villi**, which
give it the surface area it needs to absorb nutrients. In celiac disease the immune attack
flattens them. Less surface means less absorption — malnutrition, pain, fatigue, and the rest.

It is an **autoimmune disease**: the immune system attacking the body's own tissue. What makes
celiac unusual among autoimmune diseases is that the trigger is known and external. Remove
gluten and it stops.

So the real question is *why gluten* — why does wheat protein do this when rice protein does
not? The answer is a chain of six links, and the whole project sits on top of it.

### Link 1 — Grains store food in proteins

A seed needs food for its seedling, and grains store that food as **storage proteins**. Every
grain has its own family:

| Grain | Storage protein family |
|---|---|
| Wheat | gluten (gliadins + glutenins) |
| Barley | hordeins |
| Rye | secalins |
| Oat | avenins |
| Rice, corn, sorghum, millets | different ones again |

The whole group is called **prolamins**. Every grain has them; they are simply chemically
different from one another. (`data/reference/grains.csv` carries this column for all 16 grains
in the panel.)

### Link 2 — Gluten does not digest properly

When you eat protein, digestive enzymes chop it into individual **amino acids** — the 20
building blocks proteins are made of, each written as a single letter (`P`, `Q`, `E`, `F`, …).

Gluten is unusually rich in two of them: **proline (`P`)** and **glutamine (`Q`)**.

Proline is the problem. It has a rigid kink in its structure, and human digestive enzymes
cannot cut next to it. Gluten contains so much proline that the enzymes cannot finish the job.

So instead of being reduced to single amino acids, **fragments survive** — chains of roughly
10 to 33 amino acids, intact, sitting in the gut. A chain that short is called a **peptide**
(a peptide is just a short protein, under roughly 50 amino acids).

> **First coincidence:** gluten's chemistry lets fragments survive digestion at all.

### Link 3 — A gut enzyme modifies those fragments

In the gut wall there is an enzyme called **tissue transglutaminase**, abbreviated **TG2**.

TG2 takes some of the glutamines (`Q`) in those surviving peptides and converts them into
**glutamate (`E`)**. The chemical change is small; the consequence is not:

- **`Q` (glutamine)** — electrically **neutral**
- **`E` (glutamate)** — **negatively charged**

This conversion is called **deamidation**, and it is the hinge the whole disease turns on.

TG2 only acts in certain contexts: it needs a glutamine, then any amino acid, then a proline
(written **`Q-X-P`**). Gluten is loaded with both `Q` and `P`, so it carries an unusual number
of these sites.

> **Second coincidence:** gluten is an unusually good substrate for TG2.

This is why `data/reference/peptides.csv` has a `form` column holding `native` or
`deamidated`, and why 13 of the 17 panel rows are the deamidated form. The native rows (`N01`
–`N03`) exist as paired controls: same gene locus, differing only by `Q`→`E`, so the pair
isolates the effect of deamidation by itself.

### Link 4 — The immune system displays fragments for inspection

The immune system has cells that constantly pick up fragments of whatever is around and **hold
them up for inspection**.

The structure they hold fragments up *in* is an **MHC molecule** — in humans usually called an
**HLA** molecule. (MHC and HLA mean the same thing here; HLA is just the human naming system.
Do not let the two names confuse you.)

Picture a display stand with a narrow slot in it. The slot is the **binding groove**. A peptide
lies down inside it, with the ends sticking out past both openings.

A **T cell** then inspects what is on display. If it judges the fragment foreign and dangerous,
it raises an alarm and an immune response begins.

This is normal and healthy, happening thousands of times a second right now. It is how you
fight infections.

### Link 5 — Not everyone's groove is the same shape

Here is where genetics enters.

Everyone has MHC molecules, but there are many **genetic variants**, and they differ in the
shape and chemistry of the groove. Different grooves grip different peptides.

The variant that matters here is **HLA-DQ2.5** — formally the pair of genes
**DQA1\*05:01 + DQB1\*02:01**, which together build one groove.

- Roughly 25–30% of people of European descent carry it.
- **Over 90% of celiac patients carry it**; most of the rest carry a cousin, **DQ8**.
- Carrying neither makes celiac disease essentially impossible.

So DQ2.5 is necessary but not sufficient — most carriers never get sick. It is the genetic door
that has to be unlocked first.

### Link 6 — Why DQ2.5 specifically grips gluten

The floor of the binding groove is not smooth. It has small cavities along it — **pockets**.
Each one grips the side chain of whichever amino acid lands above it.

Pockets are numbered by position along the 9-amino-acid stretch that sits in the groove:
**P1, P4, P6, P7, P9**. These are the **anchor positions**, because they do most of the work of
holding the peptide in place. Amino acids at the other positions mostly point upward, out of
the groove, toward the inspecting T cell.

**DQ2.5's pockets at P4, P6 and P7 are lined with positive charge.** Positive attracts
negative, so those pockets want a negatively charged amino acid sitting in them.

Recall Link 3: TG2 has just converted neutral `Q` into negatively charged `E`.

> **Third coincidence:** TG2 installs exactly the charges DQ2.5's pockets are looking for.

This shows up directly in the panel. Every deamidated DQ2.5 core in `peptides.csv` carries a
glutamate at position 4, 6 or 7 — ten out of ten:

```
E01 PFPQPELPY  E@6     E06 PQPEQPFPQ  E@4
E02 PQPELPYPQ  E@4     E07 PFSEQEQPV  E@4,6
E03 PFPQPEQPF  E@6     E08 PQQSFPEQQ  E@7
E04 PQPEQPFPW  E@4     E11 PYPEQEEPF  E@4,6,7
E05 PIPEQPQPY  E@4     E12 PYPEQEQPF  E@4,6
```

That is not a curation artifact; it is the pocket chemistry above showing up in independently
published epitopes. It is also a cheap sanity check on any candidate peptide: a proposed DQ2.5
epitope with no route to a glutamate at an anchor position is suspect.

### The chain, assembled

> You eat wheat. Gluten's proline content means fragments survive digestion. TG2 converts some
> neutral glutamines into negatively charged glutamates. Those negative charges drop into
> DQ2.5's positively charged pockets. The peptide binds tightly, gets displayed, T cells
> attack — and the intestinal lining is damaged in the process.

Rice protein does none of this. It lacks gluten's proline content, so it digests fully; and
what survives does not have the right residues in the right places for DQ2.5's pockets.

### Framing this project holds, and keeps having to re-hold

- There is **no "gluten receptor."** The target is the **HLA-DQ2.5 binding groove**.
- You model **short peptides (9–17 aa)**, not whole proteins — peptides are what sits in the
  groove.
- Genome models are the wrong instrument. This is not a question about genomes; it is one short
  peptide fitting one slot.

---

## 2. The question this project asks

What is known clinically, from decades of feeding people things and observing the result:

- **Toxic:** wheat, barley, rye
- **Safe:** rice, corn, sorghum, millets, teff, buckwheat, quinoa, amaranth
- **Contested:** oat — most patients tolerate it, some react

That is an *observation*. Nobody derived it; it was measured in patients.

**The question is whether that pattern can be explained from the molecules alone:**

> Take peptides from toxic grains and peptides from safe grains. For each one, build a 3D model
> of it sitting in the HLA-DQ2.5 groove. Measure how well each binds. **Do the toxic ones bind
> better?**

### Why it is a good question

**It can fail.** Maybe binding strength does not separate toxic from safe. Maybe what matters
more is how *much* storage protein a grain contains, or how it is processed, or an immune
pathway nobody is modeling here. "Binding does not explain it" is a real finding, not a failed
project. A question that can only come out one way is not research.

**The answer is checkable.** The clinical pattern is already known, so there is a correct answer
to grade against — rare and valuable in computational work.

**It is the right size.** Standard tools, CPU-only, finishable on this timeline.

---

## 3. The pipeline — what each piece is for

Four things are needed for every peptide.

### 3.1 A sequence — which amino acids, in what order

Two sources:

- **`data/reference/peptides.csv`** — the panel of 17. These are **epitopes**: the specific
  short pieces the immune system is known to recognize, experimentally confirmed in the
  published literature. Each row carries a 9-amino-acid `core_9mer` (what sits in the groove)
  and a longer `extended_sequence` (the core plus real flanking residues from the source gene,
  which is what actually gets modeled).
- **`src/grains/proteomes.py`** — downloads storage proteins from all 16 grains, so **control
  peptides** can be drawn from the safe ones. These are essential. "Gluten peptides bind DQ2.5"
  means nothing without "and rice peptides do not."

### 3.2 A 3D structure — that peptide sitting in the groove

A sequence is a string of letters. Binding is physical, and depends on 3D shape. So the shape
has to be built.

**PANDORA** does this: give it a peptide sequence and an MHC allele, and it builds a predicted
3D structure of the two together. It works by **homology modeling** — finding an experimentally
solved structure of a similar complex and threading the new sequence onto that scaffold. It
uses **MODELLER** underneath, which is why a (free, academic) MODELLER license key was needed.

### 3.3 A score — how good is this fit

**PRODIGY** reads the 3D structure and counts the actual physical contacts across the interface
between peptide and groove: which atoms touch which, how many, of what kind. It returns a
number estimating binding strength.

Also computed per complex: buried surface area, anchor-pocket occupancy, and interface
confidence metrics. One row per peptide in `data/processed/scores.csv`.

### 3.4 A statistical comparison

Line the toxic-grain numbers up against the safe-grain numbers. Is the difference real or is it
noise? That is Checkpoint 4.

---

## 4. Why build structures instead of just comparing sequences

This is a real fork in the road, and worth being able to defend.

You *could* just ask "does this new peptide look like a known epitope?" That is string matching
against the catalog, and it works — for peptides similar to ones already catalogued.

But it can only ever tell you **"this resembles things already known to bind."** Its only
standard of truth is the existing list, which makes it quietly circular: the catalog is both the
teacher and the grader. It cannot speak to a peptide nobody has tested.

A structural model never consults the catalog. **PANDORA** builds geometry from physics and MHC
templates; **PRODIGY** counts atomic contacts. Neither asks whether a peptide is famous. So a
high score from this pipeline is **independent evidence** — which is the whole point when the
subject is safe grains, where no catalog of epitopes exists to match against.

> A sequence model ranks peptides by resemblance to known binders.
> A structural model ranks them by predicted binding.

### A worked example of the difference

A 2026 study (ESM-2 protein-language-model embeddings + logistic regression, per-residue
epitope prediction across grain prolamins) reported leave-one-grain-out AUC of 0.82–0.94 and
flagged a hotspot of "novel" candidates in oat avenin-E. Its own results illustrate the limit:

- **A trivial 9-mer sequence baseline matched it** (LOGO mean AUC 0.900 vs 0.903). Because the
  labels were defined by exact string match to the epitope catalog, the task *was* string
  matching — so a lookup table is the ideal tool and the language model adds nothing.
- **Its cross-species test leaked.** 74–100% of held-out positives shared a 9-mer verbatim with
  the training set, so the "unseen species" was not unseen.
- **Its top candidates land on an epitope it trained on.** Aligned against `Q09114`
  (`AVEE_AVESA`, avenin-E — the same protein panel row `E11` is drawn from), its three
  highest-scoring 9-mers all hit residues 13–19. `DQ2.5-ave-1a` (= `E11`) sits at 14–22. They
  are cultivar variants of the epitope itself, spelled with `Q` where the reference has `E`, so
  they failed to exact-match the catalog, were labeled negative, scored 0.9999, and were
  reported as discoveries.
- **Its flagship untested candidate fails the anchor check.** `PQQQPFYQQ`'s only `Q-X-P` site
  puts the glutamate at position 3 — not an anchor pocket (see §1, Link 6). *Caveat:* the
  binding **register** can shift, so this is a hypothesis rather than a verdict — and
  resolving it is exactly what a 3D model does and a sequence classifier structurally cannot.

Two things worth taking from it regardless:

- **Non-independence of repetitive prolamin sequences is a general hazard,** not a quirk of
  this panel. Its cultivar comparison collapsed because several avenin-E sequences were
  byte-identical; the equivalent here is `E03`/`E06` sharing 7 of 9 residues. See the
  Checkpoint 3 carry-forward.
- **Resample at the family level, not the peptide level.** Its protein-level block bootstrap is
  the right instinct: resampling individual positions (or peptides) pretends neighbors are
  independent and yields falsely narrow confidence intervals.

And one mistake to avoid: it set `p = 1 - classifier_score` and fed that to a
Benjamini–Hochberg correction. Classifier scores are not uniform under the null, so the
false-discovery-rate guarantee does not hold — and BH still returns a confident-looking number
anyway. For Checkpoint 4, build a real null: permute the toxic/safe labels, or score decoys
matched for **length and P/Q composition** (prolamins are compositionally extreme, so generic
decoys would flatter the result).

---

## 5. Where the project stands

Full detail in [`checkpoints.md`](checkpoints.md). The short version:

### Checkpoint 1 — get the inputs right ✅

Build the reference data; verify every row against primary literature.

**Why it mattered:** if a peptide sequence is wrong, every structure and score downstream is
wrong too — and nothing about the final number will *look* wrong. Four classes of silent error
were found and fixed: five wrong PMIDs, one row carrying a native sequence under a `deamidated`
label, five rows with fabricated flanking residues, and a readiness check that passed on an
unvalidated license.

### Checkpoint 2 — prove the pipeline reproduces a known answer ✅

Before trusting a model on unknown peptides, show it gets a *known* one right.

**1S9V** is an experimentally determined structure of DQ2.5 with a gluten peptide in the
groove, solved by **X-ray crystallography** (crystallize the protein, shine X-rays through it,
reconstruct atomic positions from the scattering pattern). That is ground truth — and
conveniently its peptide is `LQPFPQPELPY`, panel epitope `E01` exactly.

**How two structures are compared: RMSD** (root-mean-square deviation). Superpose the two,
measure how far each atom sits from where it should be, combine into one number in **ångströms**
(Å, a unit of length; a chemical bond is about 1.5 Å). RMSD 0 means identical; larger means more
different.

**The methodological core of this checkpoint:** an RMSD is meaningless on its own. Is 0.7 Å
good? There is no answer without a scale. So the scale was measured **first, from crystal
structures alone, before any model existed**:

| Comparison | E01 core RMSD | What it means |
|---|---|---|
| Two copies of the peptide in the *same* 1S9V crystal | **0.24 Å** | Noise floor — disagreement within a single experiment |
| Two *independent* crystals (1S9V vs 4OZI) | **0.91–1.16 Å** | How well the best experimental work agrees with itself |
| **AlphaFold model vs the 1S9V crystal** | **0.59–0.73 Å** | — |

Read the third row against the second: **the model is closer to the crystal than two
independent experiments are to each other.** That is the result, and it is only a result
*because* the scale was built first.

**Circularity was tested, not assumed.** 1S9V was not among AlphaFold's templates; 4OZI was. A
model merely recalling its template would hug 4OZI and sit at 4OZI's own distance from 1S9V
(~1.1 Å). Measured, it is roughly equidistant from both (0.59–0.73 Å vs 0.48–0.73 Å) while the
two crystals are 0.91–1.16 Å apart — a consensus prediction, not template recall.

**Per-residue deviation matches the structural expectation:** lowest in the 9-mer core
(0.13–0.40 Å) and highest at the flanking termini (0.83, 1.18 Å). The core is clamped in the
groove while the flanks extend into solvent, so the model is most accurate exactly where the
binding happens. The reverse pattern would have been a warning.

Figure: `results/figures/af_vs_1S9V_overlay.png`.

### Checkpoint 3 — run the real thing ⬅️ **current**, due 2026-10-20

Build structures for roughly 40–50 complexes (the 17-row panel plus matched controls), score
them all, write `data/processed/scores.csv`.

**One licensing rule to hold onto.** AlphaFold Server's terms of service forbid feeding its
output into binding- or affinity-prediction systems. PRODIGY is exactly that. So AlphaFold was
fine for Checkpoint 2 — pure geometry, no scoring — but **cannot** be used here. PANDORA is the
permitted route for anything that gets scored.

`grains.modeling.assert_scorable()` enforces this in code, default-deny, so the forbidden
combination fails loudly rather than producing a PRODIGY number that looks perfectly valid.

### Checkpoint 4 — the statistics, due 2026-10-30

Do toxic peptides score better than safe ones, and is the difference real? Plus a sequence
baseline (NetMHCIIpan via IEDB) to ask where structure adds signal beyond sequence.

### Checkpoint 5 — complete draft, due 2026-11-07

Notebook running top to bottom, `report.md`, README with a headline figure.

### 2026-11-07 → 11-28 — finishing touches and fixes

Revision, better figures, stretch goals. The draft exists from Nov 7, so this window improves
it rather than producing it.

---

## 6. The recurring failure mode

Worth stating in its own section, because it has now been paid for four times:

> **The silently plausible wrong number** — a result that arrives looking computed rather than
> erroneous.

Instances so far:

- A disordered 7-residue loop gap in 1S9V chain B that a naive CA-atom walk splices across.
- 4OZI's peptide chains being the **reverse** of their chain lettering (chain J binds A/B).
  Pairing by letter gives a 62 Å peptide RMSD — physically impossible, which is the only reason
  it was caught.
- Wrong PMIDs that still resolve to real papers.
- A readiness check that passed on an unvalidated MODELLER license.

The habit that catches them: **assert the assumption in code, and keep the bad case as a
regression test.** 50 tests, green.

---

## Glossary

| Term | Meaning |
|---|---|
| **amino acid** | One of 20 building blocks of proteins; each has a one-letter code |
| **ångström (Å)** | Unit of length, 10⁻¹⁰ m. A chemical bond is ~1.5 Å |
| **anchor position** | A peptide position (P1/P4/P6/P7/P9) whose side chain sits in a groove pocket |
| **autoimmune disease** | The immune system attacking the body's own tissue |
| **binding groove** | The slot in an MHC molecule that a peptide lies down in |
| **core (9-mer)** | The 9 amino acids of a peptide that actually sit in the groove |
| **deamidation** | TG2 converting glutamine (`Q`, neutral) to glutamate (`E`, negative) |
| **DQ2.5** | The HLA variant carried by >90% of celiac patients; DQA1\*05:01 + DQB1\*02:01 |
| **epitope** | The specific short fragment an immune receptor recognizes |
| **gliadin / glutenin** | The two protein families making up wheat gluten |
| **HLA** | Human naming system for MHC molecules; same thing |
| **homology modeling** | Building a structure by threading a sequence onto a solved similar structure |
| **IEDB** | Immune Epitope Database — community registry of observed epitopes |
| **MHC** | Major histocompatibility complex; the display molecule holding peptides up for T cells |
| **MODELLER** | Homology-modeling program; PANDORA's engine. Needs a free academic license |
| **PANDORA** | Builds peptide–MHC 3D models by homology. CPU-only |
| **peptide** | A short protein, under roughly 50 amino acids |
| **pocket** | A cavity in the groove floor that grips one amino acid side chain |
| **PRODIGY** | Scores a structure's interface to estimate binding strength |
| **prolamin** | The seed storage protein group; gluten, hordein, secalin, avenin, zein, … |
| **register** | Which 9 residues of a longer peptide actually occupy the groove |
| **RMSD** | Root-mean-square deviation; one number for how different two structures are |
| **T cell** | Immune cell that inspects what MHC molecules display |
| **TG2** | Tissue transglutaminase; the gut enzyme that deamidates gluten peptides |
| **villi** | Finger-like projections lining the small intestine; flattened in celiac disease |
| **X-ray crystallography** | Experimental method for determining atomic structure |
