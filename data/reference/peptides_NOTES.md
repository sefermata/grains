# Peptide panel — verification notes

`peptides.csv` is a **working draft**. Before Checkpoint 1 is "done", every row's `verified`
field must read `yes`, cross-checked against a primary source.

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

Still open for Checkpoint 1: add `iedb_id` column (checklist step 2), add UniProt accession +
residue range column (checklist step 3), resolve the E10 name-suffix question, and ideally
find an accessible copy of Tye-Din 2010 to table-check N01–N03 directly rather than by
derivation.

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
