# Manual steps

Two parts of the pipeline can't be automated: MODELLER needs a license key tied to a person,
and AlphaFold Server has no public API. Both are one-time (or near-one-time) and both block
Checkpoint 2/3, so the exact steps live here rather than in a chat log.

---

## 1. MODELLER + PANDORA  (blocks Checkpoint 3 — the critical path)

PANDORA is the **only** modeling route allowed to feed the scoring step (AlphaFold Server's
Prohibited Use Policy bars its output from binding predictors — see `src/grains/modeling.py`).
PANDORA is a thin layer over MODELLER, and MODELLER needs a free academic license key.

### Step 1 — register for the key

<https://salilab.org/modeller/registration.html>

- Platform: check **Apple Mac OS X** only. Nothing else applies, and *don't* check the
  web-service option (ModWeb/Chimera) — we run MODELLER locally.
- The platform choice only decides which download links you're shown. Since the install comes
  from conda, the download is irrelevant; **the key is the point**, and it's platform-independent.
- Email: use an **institutional address** (e.g. `@sierracollege.edu`) if you have one. The
  academic license asks for one and free webmail is often rejected.

The key arrives by email — a short word-like string.

### Step 2 — install

```bash
conda install -y -c salilab modeller      # 10.8 has a py313 build; ~200 MB
pip install csb-pandora                   # 0.9 on PyPI; imports MODELLER, so install it second
```

### Step 3 — enter the key

conda's installer may prompt for it. If it doesn't, edit the key into MODELLER's config:

```bash
# find the file (don't guess the path — the version number is in it)
find /opt/anaconda3 -name config.py -path '*modeller*'
# then replace the XXXX placeholder on the `license` line with your key
```

### Step 4 — verify

```bash
cd ~/grains
PYTHONPATH=src python -c "from grains.modeling import pandora_available; print(pandora_available())"
```

Expect `(True, 'PANDORA and MODELLER available')`. Anything else prints the reason.

### Watch the disk

This Mac is at ~94% (13 GB free). MODELLER itself is modest, but PANDORA builds a **local
template database** of MHC structures on first use, which is the part that can get large.
Check free space before and after:

```bash
df -h /Users/sefoster | tail -1
```

If the database build threatens the remaining space, stop and reconsider scope — a smaller
template set is a legitimate choice to make and state, and it's cheaper than filling the disk.

---

## 2. AlphaFold Server  (blocks the Checkpoint 2 RMSD comparison)

Used for **exactly one** structure: the validation complex. Its only permitted downstream use
is a geometric RMSD against the 1S9V crystal. **Never a PRODIGY input** — that's the ToS
constraint, enforced in code by `modeling.assert_scorable()`.

### Step 1 — generate the submission file (already done, regenerate any time)

```bash
cd ~/grains
PYTHONPATH=src python -m grains.modeling
# writes data/processed/alphafold_validation_job.json
```

Three protein chains: DQA1\*05:01 (193 aa), DQB1\*02:01 (198 aa), peptide `LQPFPQPELPY`.
The peptide is the sequence as deposited in 1S9V chain C, so the model and the crystal are the
same molecule and the RMSD needs no sequence caveat. No PTM is declared for the deamidated
residue: TG2 converts Gln to Glu, Glu is a standard amino acid, so it's just the `E`.

### Step 2 — submit

1. <https://alphafoldserver.com> — sign in with a Google account.
2. Import the JSON (the bulk-upload / "import JSON" control). It loads as a **draft** job.
3. Run it. The free tier allows 20 jobs/day, so one validation job costs nothing scarce.

### Step 3 — what to bring back

Download the result and unzip it into `data/raw/models/`:

```bash
mkdir -p data/raw/models
# unzip the downloaded archive there
```

**The file that matters is the predicted structure — the `.cif` (mmCIF).** That's the only
thing the RMSD needs. The confidence JSON (pLDDT/PAE) is worth keeping for the writeup as a
record of how confident the prediction was, but note that interface PAE and peptide pLDDT in
`docs/checkpoints.md`'s Checkpoint 3 scoring list apply to the **PANDORA** models, not this one.

Then tell me it's there, and I'll superpose it against both 1S9V copies and report the RMSD
next to the crystal-vs-crystal baseline (E01 core: 0.24 Å within one crystal, ~1.2 Å between
independent depositions — that's the bar).

`data/raw/` is git-ignored, so the downloaded model stays local.
