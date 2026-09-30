# Tutor Center Project — Rationale and Alignment

> Companion to the project itself. Addresses the Tutor Center project rubric:
> **Audience**, **Purpose and Relevance** (including CRLA alignment, high-impact practices,
> and sources). Items marked **[confirm]** are Sierra College / Tutor Center specifics to
> verify against current information before submission — they are not asserted in any
> student-facing material.

## The project in one paragraph

A five-part, self-contained data-science project, *Spotting Gluten Grains from Protein
Sequence*. Students use only introductory data-science skills — tables, grouped
visualization, a permutation test, and a k-nearest-neighbors classifier — to
test whether the amino-acid makeup of a grain's storage proteins can separate the grains
that cause celiac disease (wheat, barley, rye) from those that do not. It ships as paired
Jupyter notebooks (a Berkeley `datascience` track and a `pandas` track), a full solution
set, an instructor guide, a dataset that regenerates from a script, and a held-out "mystery
set" students classify at the end. Because the toxic/safe status of every grain is already
known, the project has a built-in answer key.

---

## Audience

**Primary audience.** Students in the Sierra College introductory data-science course
([course number and title — confirm], modeled on UC Berkeley's *Foundations of Data
Science* / "Data 8") who come to the Tutor Center for help with that course. These students
are typically:

- new to programming, with exposure limited to the course's own Python and `datascience`
  environment;
- not biology majors — no biology or chemistry coursework is assumed;
- most engaged by the course's prediction unit, and by tasks that have a concrete,
  checkable answer rather than an open-ended write-up.

**Secondary audience.** Tutor Center tutors who support this course and want a complete,
worked example to reference during the classification unit; Tutor Center coordinators
reviewing discipline-specific materials for the project collection.

**How the audience shaped the design.** (This is the "Style and Tone" consideration.)

| Audience need | Design response |
|---|---|
| No biology background | A short "The biology you need (all of it)" section; every term — protein, amino acid, proline, `accession`, feature, baseline — defined at first use |
| Still building Python fluency | Scaffolded code cells (function signature plus blanks), not empty cells; the canonical solution to the counting question is an explicit `for` loop, with the one-line version only mentioned |
| Drop-in setting; a tutor may not be beside the student | A `check_*` helper after **every** question prints `OK` or a specific hint — immediate formative feedback |
| Course uses `datascience` now, but students move on to `pandas` | Two parallel notebooks; the `pandas` track adds "pandas note" callouts at each point the two libraries diverge |
| Reflection without overload | Written prompts are one sentence, with a concrete computational lead-in (compute the numbers first, then interpret) |
| American English convention | American spelling and plain second-person voice throughout |

---

## Purpose and Relevance

**The need it addresses.** The course's built-in projects are well constructed but generic
(movie ratings, climate data). Two gaps remain: students rarely see where these skills lead
in a STEM research context, and tutors have no ready end-to-end worked example for the
classification unit. This project fills both. It is:

- **authentic** — a real, published research question in miniature (structural immunology of
  gluten), reduced to something an intro student can actually do;
- **self-contained** — no licensed software, no internet needed at runtime; the dataset
  rebuilds from `build_dataset.py` and the notebooks from `build_notebooks.py`;
- **self-checking** — the known toxic/safe grain list is the answer key, so a student can
  judge their own result, and the "mystery set" gives an honest held-out test.

**Use in the current semester.**

- Offered at the Tutor Center as an optional enrichment project during the
  prediction / classification unit ([course weeks — confirm]).
- Run as a one-hour tutor-led group session: work Parts 0–2 together, assign Parts 3–5.
- Used at tutor training ([session — confirm]) to demonstrate questioning strategy and
  modeling of problem-solving on Part 5.

**Use in future semesters.** Nothing is tied to a term. The biology is stable; the build
scripts regenerate every artifact; the `datascience` + `pandas` split protects against a
change in the course's tooling. The project can live in the Tutor Center project collection
and be reused or adapted (e.g., a shorter Parts 0–3 lab).

**Where / when / how / to whom.** Tutor Center drop-in and appointment tutoring, and group
tutoring sessions, during the classification unit; presented to the tutor cohort at
training; archived for future tutors.

**Alignment with CRLA training.** Preparing this project, and using it with tutees, exercises
several College Reading & Learning Association (CRLA) International Tutor Training Program
topic areas:

- **Questioning techniques / Socratic method** — Q0.3 ("find the score to beat") and Q1.4
  ("which number is a fair comparison, and why?") are written as guided questions, not
  explanations.
- **Modeling problem-solving** — the solution notebooks lay out explicit step-by-step
  reasoning rather than terse idioms.
- **Active and inductive learning** — students discover the signal (gluten proteins are
  proline/glutamine-rich) and, later, its limits (maize and sorghum proteins look
  suspiciously similar) from the data themselves.
- **Fostering independent learners / metacognition** — the `check_*` feedback loop and the
  one-sentence reflections build self-monitoring.
- **Transfer of learning** — the two-track design explicitly bridges the classroom tool
  (`datascience`) and a professional one (`pandas`, and `scikit-learn` in Part 5).

**High-impact practices.** For an introductory student, the project approximates an
undergraduate-research experience — an authentic question, real data from a public
repository, an honest evaluation on held-out data, and a written conclusion — and, run as a
group session, supports collaborative learning.

**Academic, audience-centered content.** Every biological claim in the student materials is
traceable to the peer-reviewed literature listed below. The dataset is derived from UniProt,
the standard public protein-sequence repository. No claim about Sierra College resources
appears in any student-facing file; the administrative details in *this* document are marked
**[confirm]** and must be checked against current Tutor Center and course information.

---

## Sources

APA 7th edition. **Verify each DOI / PMID against the publisher before final submission.**

**Primary literature (basis for the biology):**

- Shan, L., Molberg, Ø., Parrot, I., Hausch, F., Filiz, F., Gray, G. M., Sollid, L. M., &
  Khosla, C. (2002). Structural basis for gluten intolerance in celiac sprue. *Science,
  297*(5590), 2275–2279. https://doi.org/10.1126/science.1074129
- Tye-Din, J. A., Stewart, J. A., Dromey, J. A., Beissbarth, T., van Heel, D. A., Tatham,
  A., Henderson, K., Mannering, S. I., Gianfrani, C., Jewell, D. P., Hill, A. V. S.,
  McCluskey, J., Rossjohn, J., & Anderson, R. P. (2010). Comprehensive, quantitative
  mapping of T cell epitopes in gluten in celiac disease. *Science Translational Medicine,
  2*(41), 41ra51. https://doi.org/10.1126/scitranslmed.3001012
- Sollid, L. M., Qiao, S.-W., Anderson, R. P., Gianfrani, C., & Koning, F. (2012).
  Nomenclature and listing of celiac disease relevant gluten T-cell epitopes restricted by
  HLA-DQ molecules. *Immunogenetics, 64*(6), 455–460. PMID 22371225.
  https://doi.org/10.1007/s00251-012-0599-z
- Sollid, L. M., Tye-Din, J. A., Qiao, S.-W., Anderson, R. P., Gianfrani, C., & Koning, F.
  (2020). Update 2020: Nomenclature and listing of celiac disease-relevant gluten epitopes
  recognized by CD4+ T cells. *Immunogenetics, 72*(1–2), 85–88.
  https://doi.org/10.1007/s00251-019-01141-w

**Data source:**

- The UniProt Consortium. (2023). UniProt: The Universal Protein Knowledgebase in 2023.
  *Nucleic Acids Research, 51*(D1), D523–D531. https://doi.org/10.1093/nar/gkac1052

**Tools and course environment:**

- Adhikari, A., DeNero, J., & Wagner, D. (2022). *Computational and inferential thinking:
  The foundations of data science* (2nd ed.). https://inferentialthinking.com
- The pandas development team. (2024). *pandas documentation.*
  https://pandas.pydata.org/docs/
- Pedregosa, F., Varoquaux, G., Gramfort, A., Michel, V., Thirion, B., Grisel, O., Blondel,
  M., Prettenhofer, P., Weiss, R., Dubourg, V., Vanderplas, J., Passos, A., Cournapeau, D.,
  Brucher, M., Perrot, M., & Duchesnay, É. (2011). Scikit-learn: Machine learning in
  Python. *Journal of Machine Learning Research, 12*, 2825–2830.

The career framing in the notebook's "Why this project" section is motivational and
descriptive (data-analysis roles in biology have grown alongside the size of biological
datasets); it does not cite a specific labor statistic and should not be read as one.
