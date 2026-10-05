# Peptide panel — verification notes

`peptides.csv` is a **working draft**. Before Checkpoint 1 is "done", every row's `verified`
field must read `yes`, cross-checked against a primary source.

## IEDB lookup, 2026-10-05 — checklist step 2 done

`grains.iedb` resolves each core to an IEDB epitope id; `iedb_id` is now a column in
`peptides.csv`, with the full result (all candidates per row) in
`data/processed/peptide_iedb_map.csv`. **12/15 cored rows resolved.** Responses are cached on
disk, so re-runs cost no requests.

Two lookups, because neither covers the panel alone. An **exact match on `core_9mer`** hits
when the epitope is registered as the bare 9-mer (E01 → 47539). Otherwise a **contains-search**
finds the longer peptide that was actually synthesised and assayed, with our core inside it.
Exact search on `extended_sequence` finds nothing for any row — our window is cut from the gene,
not copied from whatever length a given lab used.

**Candidate choice is specificity-first:** shortest registered peptide containing the core,
then most assays, then lowest id. The first attempt ranked by assay count instead and gave E03
*and* E06 the same id — IEDB 462152, a 24-residue ω-gliadin repeat
(`PQPEQPFPQPEQPFPQPEQPFPQP`) containing both cores. An id shared between two rows identifies
neither, so the rule now prefers the most specific entry. A duplicate-id check runs regardless
and writes a warning into the output.

### E03 and E06 are nearly the same peptide — a real result, not a lookup artifact

The duplicate warning still fires for E03/E06, and investigating it showed the cause is
biological:

- **DQ2.5-glia-ω1** (E03, wheat ω-gliadin) core `PFPQPEQPF`
- **DQ2.5-sec-2** (E06, rye secalin) core `PQPEQPFPQ`

These **overlap by 7 of 9 residues** (`PQPEQPF`) — two offset registers of the same repetitive
prolamin motif. Every one of E06's three IEDB candidates also contains E03's core; no
registered peptide separates them. So the shared id is correct and unavoidable, and it is
recorded with a warning rather than papered over.

Worth noting this is part of *why rye is toxic*: rye secalin carries essentially the same
DQ2.5 epitope motif as wheat ω-gliadin. That is a point for the writeup, not a problem.

**Carry forward to Checkpoints 3–4:** E03 and E06 will model and score almost identically, so
they are **not independent observations**. Counting both toward a toxic-vs-safe comparison
overstates the sample. Either treat them as one observation, or report the statistics with and
without E06 and say which was used.

**Three rows unresolved** — E07 (DQ2.5-glut-L1, LMW-glutenin), E11 and E12 (the two avenins).
Their cores are not registered in IEDB alone or inside any longer peptide. These are the less
studied members of the panel, so absence from the database is plausible rather than suspicious;
left blank and stated rather than filled with a near-miss.

## UniProt mapping, 2026-10-05 — checklist steps 3 and 4 done

`grains.map_peptides` locates every peptide in the downloaded proteomes and writes
`data/processed/peptide_uniprot_map.csv`. `uniprot_accession` and `uniprot_range` columns are
now in `peptides.csv`. **17/17 rows are gene-faithful**: 13 deamidated rows each need at least
one Q→E against the gene, 4 native rows match literally. That consistency is itself the check.

**Step 3 as written doesn't survive the data.** A deamidated peptide is *not* a substring of
any native gene — that is what deamidation means. Searching literally reports every deamidated
epitope as absent. So the search lets a peptide `E` match a gene `Q` (never the reverse: TG2
goes Q→E), one degree of freedom per glutamate with every other residue literal, and records
which positions flexed. That does step 3 and step 4 at once: a position that had to flex from
Q *is* a deamidation site.

**What it confirmed.** E01/N01 land on the identical window of P04724 (77–93) differing by
exactly one Q→E; likewise E02/N02 and E04/N03; and C01 flexes at exactly 3 positions against a
literal C02. The native/deamidated pairs are controlled pairs at one locus — established from
sequence, independently of the literature.

**Two errors it caught**, both of the silent kind:

1. **E09 carried a native `extended_sequence`** (`QGSFQPSQQNPQ`) while labeled
   `form=deamidated`, and it didn't even contain its own deamidated core `EGSFQPSQE`. It was
   the only deamidated row matching a gene literally — zero Q→E — which is what gave it away.
   Since `extended_sequence` is the modeling input, E09 would have been modeled native under a
   deamidated label, inverting the TG2 comparison it exists to support. Corrected to
   `EGSFQPSQENPQ`, the C-terminal 12 residues of 2NNA's deposited peptide.
2. **Five rows had fabricated flanks** (E05, E06, E08, E11, E12). Each verified core sits in a
   real gene, but the surrounding residues of `extended_sequence` matched no gene at all —
   residues appear to have been dropped at a junction (E05's gene reads `FPQQ|PIPQQPQPY`, but
   the row spliced `F` straight onto `PIPEQ`). E12 is the clearest: it ended `PFVQ`, but `V`
   belongs to E11's paralog; the gene holding E12's core has `PFLQQ`. All five replaced with
   real gene windows.

**First diagnosis of those five was wrong, and worth recording as a lesson.** The initial run
blamed the rye/teff/oat proteome coverage gap. That was the convenient explanation and it was
false: the cores *are* present locally, including in rye (`A4GU92_SECCE`) and oat
(`AVEE_AVESA`, `F4MJY2_AVESA`). The mapper now checks for the core before blaming coverage and
reports `match_type=core_only` with the real location, so the diagnosis points at the panel
rather than at the data.

### Two curation decisions, and why

**Which paralog, when a core occurs in many?** Prefer **reviewed (Swiss-Prot)** entries over
unreviewed (TrEMBL), then lowest accession to break ties. Storage-protein families are large
and repetitive — E08's core appears in **101** wheat entries with differing flanks — so this
needs a rule, not a preference. Swiss-Prot entries are manually curated with stable names and
sequences, which makes them defensible to cite; the accession tiebreak keeps the choice
deterministic so output doesn't shuffle between runs. This moved E08 to `P06659 GDBB_WHEAT`
("Gamma-gliadin B", Swiss-Prot). Two rows have **no reviewed entry available**: E06 (rye) and
E12 (oat) rest on TrEMBL — a real limitation of public annotation for those grains, stated
rather than hidden.

**How long should `extended_sequence` be?** Each repaired row **keeps its existing length**
(13–14 aa); only the flanks changed. Rejected the tidier option of standardizing the panel to
one length, for two reasons. The rows that were already gene-faithful (E01–E04, E07, E09, E10)
carry their extended sequences as published, and rewriting those to hit a round number would
trade citable sequences for a neater table. And peptide length drives the Checkpoint 3
interface scores directly — a longer peptide buries more surface whatever its chemistry — so
what matters is that the panel's length profile is *known* and that controls are matched to it
peptide-by-peptide, which `build_controls.py` already plans to do. Repairing these five at
their current lengths keeps that profile intact instead of stacking a second change on top of
the flank fix.

**Carry forward to Checkpoint 4:** panel `extended_sequence` lengths run 12–18 aa (plus the two
33-mers, which are a separate `toxic_context` category). Length is therefore a confound for any
raw interface score. Handle it by per-peptide length-matched controls, and report length as a
covariate or use per-residue normalized scores — don't compare raw buried surface area across
peptides of different lengths.

## Citation audit, 2026-10-05 — three more wrong PMIDs fixed

Every PMID in the repo was checked against PubMed `esummary` (title + journal + authors), not
spot-checked. **Three of the eleven were wrong**, including the single most-cited one:

| Was | Actually points to | Corrected to |
|-----|--------------------|--------------|
| `22371225` | a Ti/Cu-catalysis organic-chemistry paper, *Org Biomol Chem* 2012 | **`22322673`** — Sollid 2012 |
| `20679571` | "Herpes zoster vaccine in nursing facility residents", *Ann Intern Med* 2010 | **`20650871`** — Tye-Din 2010 |
| `31735986` | "Location of colorectal adenomas and serrated polyps…", *Int J Colorectal Dis* 2019 | **`31735991`** — Sollid 2020 |

`22371225` was in `source_pmid` for **all 12 epitope rows** and `20679571` for 7 more.

**The verification work itself stands.** The 2026-09-16 pass was done against **PMC3349865**,
and PMC3349865 *is* Sollid et al. 2012, "Nomenclature and listing of celiac disease relevant
gluten T-cell epitopes restricted by HLA-DQ molecules", *Immunogenetics* — confirmed by
`esummary`, which also gives its real PMID as 22322673. So the cores were checked against the
right table; only the number recorded beside them was wrong. No epitope sequence changes.

The other eight PMIDs verified correct: 12351792 (Shan 2002), 14517794 (Vader 2003), 15020763
(Kim 2004), 15526039 (Arentz-Hansen 2004), 17629515 (Henderson 2007), 23063329 (Broughton
2012), 24777060 (Petersen 2014) — plus 12923521 and 20951824, which were already correctly
identified as the *wrong* ones in the 2026-09-16 pass.

**Lesson: PMIDs recorded from memory are unreliable and fail silently** — a wrong 8-digit
number still resolves to a real paper, so nothing looks broken. Five of eleven were wrong
across the two audits. Any new citation gets `esummary`-checked before it lands; the one-liner
is in the Sources section below.

## Status as of 2026-09-16

All 17 rows now read `verified=yes` — **step 1 of the checklist below** (canonical name,
core_9mer, restriction confirmed against Sollid 2012 Table 2, PMC3349865) is done for every
row. Steps 2–3 (IEDB epitope ID, UniProt accession + residue range) are **still open** —
`iedb_id` and UniProt-accession columns don't exist in `peptides.csv` yet. Don't read
`verified=yes` as "fully done per this checklist"; it means "name/core/restriction confirmed
against a primary source."

What changed in this pass:
- E05–E09, E11, E12 (`partial`/`no`) — core_9mer confirmed exact against Sollid 2012 Table 2
  for all of them, including E08 (gamma-gliadin) which was flagged as an unverified risk.
- E10 — core confirmed, but the table extraction rendered the name as "DQ8-glia-γ1a" (with an
  "a" suffix) vs this row's "DQ8-glia-g1". Not renamed without a closer look at the primary
  text — flagged in the row's notes.
- N01–N03 (native pairs) — confirmed by reverting the single validated TG2 site (E→Q) in the
  already-Sollid-confirmed deamidated core; **not** independently checked against the
  Tye-Din 2010 table, which has no open-access full text. Weaker verification than the others
  — algebraically consistent, not table-cross-checked.
- C01/C02 (33-mer) — sequence cross-checked against a literature-quoted native 33-mer; C01's
  triple Q→E pattern matches the expected embedded epitope cores.
- **Two citation errors caught and fixed:** C01/C02 cited PMID `12923521` for Shan et al. 2002
  — that PMID actually belongs to an unrelated orthodontics-ethics paper; corrected to
  `12351792` (confirmed via PubMed eutils). E11/E12 cited PMID `20951824` for the avenin
  source — that PMID belongs to an unrelated nitric-oxide/brain-perfusion paper; corrected to
  `15526039` (Arentz-Hansen et al. 2004, PLoS Medicine — the paper Sollid 2012 itself cites for
  these two epitopes, alongside Vader et al. 2003 Gastroenterology, PMID `14517794`).

Still open for Checkpoint 1 (updated 2026-10-05): add `iedb_id` column (checklist step 2) and
resolve the E10 name-suffix question. Checklist step 3 (UniProt accession + residue range) and
step 4 (deamidation positions) are **done** — see the UniProt mapping section above.
Tye-Din 2010 and Sollid 2020 are confirmed paywalled (no PMC record), so N01–N03 stay
derivation-verified; that is now a stated limitation rather than an open task. Note the
mapping independently confirmed N01–N03 sit at the same gene loci as their deamidated
partners, differing only by Q→E, which is stronger evidence than the original derivation.

## What each column means

- `core_9mer` — the 9-residue segment that sits in the HLA groove (P1–P9). For DQ2.5 the
  key anchor positions are P1, P4, P6, P7, P9.
- `extended_sequence` — the longer peptide actually synthesised/modeled (the groove needs
  flanking residues; ~13–17 aa is typical). This is what goes to AlphaFold Server.
- `form` — `deamidated` = after TG2 (specific Q→E); `native` = as translated. DQ2.5 strongly
  prefers deamidated; `native_pair` rows exist to measure that effect.
- `immunodominant` — flagged `yes` for the epitopes shown to dominate the T-cell response in
  Tye-Din et al. 2010 (PMID 20650871): glia-α1a, glia-α2, glia-ω1, glia-ω2.

## Verification checklist (Checkpoint 1)

1. Open Sollid LM et al., *Immunogenetics* 2012 (PMID 22322673) **and** the 2020 update
   (PMID 31735991). Confirm each `canonical_name`, `core_9mer`, and `restriction`.
2. For each epitope, pull the IEDB entry (search iedb.org by name or sequence) and record the
   IEDB epitope ID in a new `iedb_id` column.
3. Confirm the `extended_sequence` is a real substring of a UniProt entry for that protein +
   grain (the `proteomes.py` output). Record the UniProt accession + residue range.
4. Confirm the deamidation positions: which Q became E, and that a real TG2 QxP-type motif
   supports it (see `deamidation` rules to be added at Checkpoint 4).
5. ~~Rows currently `verified = no`~~ — resolved 2026-09-16, see "Status" section above.
6. ~~Rows `verified = partial`~~ — resolved 2026-09-16, see "Status" section above.

## Negative / control peptides (added at Checkpoint 3, not here)

Controls come from seed-storage proteins of **safe** grains — zein (maize), kafirin (sorghum),
orzenin/glutelin (rice), plus a poly-Pro/Gln stretch from wheat gluten that is *not* a known
epitope. They will be:

- length-matched to the toxic panel (14–17 aa),
- matched on proline + glutamine fraction (so the contrast is not just "less P/Q"),
- deamidated in silico with the same TG2 rules,
- generated by `src/grains/build_controls.py` and written to `data/processed/controls.csv`.

## Sources

**Check any new PMID before recording it** (a wrong one resolves to a real paper and so fails
silently — see the citation audit above):

```bash
curl -sS "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&id=PMID&retmode=json" \
  | python3 -c "import sys,json;d=json.load(sys.stdin)['result'];k=d['uids'][0];print(d[k]['title']);print(d[k]['fulljournalname'],d[k]['pubdate'])"
```

- PMID 22322673 — Sollid et al., gluten epitope nomenclature, *Immunogenetics* 2012 (PMC3349865,
  open access — used directly for the 2026-09-16 verification pass). **Corrected 2026-10-05**
  — was miscited as PMID 22371225 (an organic-chemistry paper) in all 12 epitope rows.
- PMID 31735991 — Sollid et al., nomenclature update, *Immunogenetics* 2020 (not yet consulted
  directly — the 2026-09-16 pass used the 2012 paper only). **Paywalled**: no PMC record.
  **Corrected 2026-10-05** — was miscited as PMID 31735986 (a colorectal-adenoma paper).
- PMID 20650871 — Tye-Din et al., "Comprehensive, quantitative mapping of T cell epitopes in
  gluten in celiac disease", *Sci Transl Med* 2010. **Paywalled**: confirmed no PMC record, so
  N01–N03 stay verified by derivation only — a real, stated limitation, not an oversight.
  **Corrected 2026-10-05** — was miscited as PMID 20679571 (a herpes-zoster-vaccine paper).
- PMID 12351792 — Shan et al., structural basis for gluten intolerance (the 33-mer), *Science*
  2002. **Corrected 2026-09-16** — was previously miscited as PMID 12923521, which belongs to
  an unrelated orthodontics-ethics article.
- PMID 15526039 — Arentz-Hansen et al., the molecular basis for oat intolerance, *PLoS Medicine*
  2004. **Corrected 2026-09-16** — was previously miscited as PMID 20951824, which belongs to
  an unrelated nitric-oxide/brain-perfusion article. Sollid 2012 cites this paper (and Vader et
  al. 2003 *Gastroenterology*, PMID 14517794) for DQ2.5-ave-1a/1b.
