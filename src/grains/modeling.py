"""Two routes from sequence to a peptide-HLA complex structure, kept strictly separate.

The project needs modeled complexes for two different jobs, and a licensing constraint means
one tool may not do both:

  ROUTE_ALPHAFOLD ("alphafold_server")
      AlphaFold Server, submitted by hand through the web UI (there is no public API).
      Used for **exactly one** structure: the Checkpoint 2 validation complex, whose only
      downstream use is a geometric RMSD comparison against the 1S9V crystal structure.

  ROUTE_PANDORA ("pandora")
      PANDORA (CPU-only, MODELLER-based homology modeling). The route for every complex
      that will be scored -- the full ~40-50 peptide panel in Checkpoint 3.

Why the split: the AlphaFold Server Prohibited Use Policy (clause 2) forbids using its
output "in connection with any automated system that predicts the binding or interaction of
the protein with ligands or peptides". PRODIGY is precisely that kind of system -- the same
category as the policy's own named examples (Glide, AutoDock). An RMSD is a pure geometric
comparison and is fine; a binding score is not. `assert_scorable()` below makes that a
runtime check rather than a rule someone has to remember, because the failure mode is
invisible: a PRODIGY number computed from an AlphaFold model looks exactly like a valid one.

If PANDORA's accuracy on the validation peptide is too poor to trust for the panel, that is
a limitation to report plainly -- not a reason to reach for AlphaFold Server here instead.

Usage
-----
    PYTHONPATH=src python -m grains.modeling

Writes the AlphaFold Server submission file for the validation complex. Upload it at
alphafoldserver.com ("Bulk upload" / import JSON), run the job, then download the
predicted structure to data/raw/models/ and compare it with:

    PYTHONPATH=src python -m grains.superpose     # for the crystal-vs-crystal baseline
"""

from __future__ import annotations

import contextlib
import io
import json
import time
from pathlib import Path

from Bio import SeqIO

from . import config

ROUTE_ALPHAFOLD = "alphafold_server"
ROUTE_PANDORA = "pandora"

# Only PANDORA models may be fed to a binding/affinity scorer. See the module docstring.
SCORABLE_ROUTES = frozenset({ROUTE_PANDORA})

HLA_FASTA = config.PROCESSED / "hla_dq25.fasta"
MODELS_DIR = config.RAW / "models"

# The validation peptide: the sequence as deposited in 1S9V chain C, so the model and the
# crystal structure are the same molecule and the RMSD needs no sequence caveat. Its core
# PFPQPELPY is panel epitope E01 (DQ2.5-glia-a1a).
#
# No post-translational modification needs to be declared to AlphaFold Server for this
# peptide even though it is the deamidated form: TG2 deamidation converts glutamine to
# glutamate, and glutamate is a standard amino acid, so the deamidated residue is simply
# the E at position 9 of the sequence below.
VALIDATION_PEPTIDE = "LQPFPQPELPY"
VALIDATION_PEPTIDE_EPITOPE = "E01"


def assert_scorable(route: str, context: str = "") -> None:
    """Raise unless `route` may legitimately feed a binding-affinity scorer.

    Call this at the top of any function that computes a binding score (PRODIGY dG/Kd and
    anything equivalent) from a modeled structure. Checkpoint 3's scoring step runs this
    per complex.
    """
    if route in SCORABLE_ROUTES:
        return
    where = f" ({context})" if context else ""
    raise ValueError(
        f"route {route!r} must not be used for binding/affinity scoring{where}. "
        "The AlphaFold Server Prohibited Use Policy (clause 2) forbids using its output "
        "with automated systems that predict protein-ligand or protein-peptide binding, "
        "which is what PRODIGY does. Model this complex with PANDORA "
        f"(route {ROUTE_PANDORA!r}) instead; AlphaFold Server output may only be used for "
        "geometric comparison, e.g. the Checkpoint 2 RMSD validation."
    )


def load_hla_chains(fasta_path: Path = HLA_FASTA) -> tuple[str, str]:
    """The DQA1*05:01 and DQB1*02:01 ectodomain sequences, in that order.

    Produced by `grains.fetch_hla` from 1S9V's SEQRES records, so these are the same
    constructs the model will be superposed against.
    """
    if not fasta_path.exists():
        raise FileNotFoundError(
            f"{fasta_path} not found -- generate it first with "
            "`PYTHONPATH=src python -m grains.fetch_hla`"
        )
    records = list(SeqIO.parse(fasta_path, "fasta"))
    if len(records) != 2:
        raise ValueError(f"expected 2 records in {fasta_path}, found {len(records)}")
    alpha, beta = (str(record.seq) for record in records)
    return alpha, beta


def alphafold_server_job(peptide: str, name: str, *, alpha: str, beta: str) -> dict:
    """One AlphaFold Server job dict: the DQ2.5 heterodimer plus the peptide, 3 chains.

    Schema per AlphaFold Server's own JSON format ("alphafoldserver" dialect, version 1):
    a list of job dicts, each with `name`, `modelSeeds`, `sequences`, `dialect`, `version`.
    An empty `modelSeeds` lets the server choose the seed.
    """
    return {
        "name": name,
        "modelSeeds": [],
        "sequences": [
            {"proteinChain": {"sequence": alpha, "count": 1}},
            {"proteinChain": {"sequence": beta, "count": 1}},
            {"proteinChain": {"sequence": peptide, "count": 1}},
        ],
        "dialect": "alphafoldserver",
        "version": 1,
    }


def write_validation_submission(out_path: Path | None = None) -> Path:
    """Write the single-job AlphaFold Server submission for the Checkpoint 2 validation."""
    config.ensure_dirs()
    alpha, beta = load_hla_chains()
    job = alphafold_server_job(
        VALIDATION_PEPTIDE,
        name=f"DQ2.5_{VALIDATION_PEPTIDE_EPITOPE}_validation_1S9V",
        alpha=alpha,
        beta=beta,
    )
    out_path = out_path or config.PROCESSED / "alphafold_validation_job.json"
    out_path.write_text(json.dumps([job], indent=2) + "\n")

    print(f"wrote {out_path}")
    print(f"  job name: {job['name']}")
    print(f"  chain 1 (DQA1*05:01): {len(alpha)} aa")
    print(f"  chain 2 (DQB1*02:01): {len(beta)} aa")
    print(f"  chain 3 (peptide {VALIDATION_PEPTIDE_EPITOPE}): "
          f"{VALIDATION_PEPTIDE} ({len(VALIDATION_PEPTIDE)} aa)")
    print()
    print("Next: upload at alphafoldserver.com, run the job, download the predicted")
    print(f"structure (mmCIF) into {MODELS_DIR}/, then superpose it on 1S9V.")
    print("Its output is for the RMSD check ONLY -- never a PRODIGY input (see docstring).")
    return out_path


def pandora_available() -> tuple[bool, str]:
    """Whether the PANDORA route can run here, and why not when it can't.

    PANDORA needs MODELLER, which is free for academic use but requires a license key
    (registration at salilab.org/modeller/registration.html).

    The license is checked by constructing a `modeller.Environ()`, not by importing
    modeller. Importing is not a sufficient test: MODELLER holds its key in
    `modlib/modeller/config.py`, and whether a bad key fails at import or only at first
    use is not something to rely on -- `Environ()` is the call that actually validates it.
    A readiness check that answers "ready" against an unusable license is worse than no
    check, because the failure then surfaces somewhere further down the pipeline.
    """
    try:
        import modeller
    except ImportError as exc:
        return False, (
            "MODELLER is not installed. Install it and PANDORA with:\n"
            "    export KEY_MODELLER='your-key'   # before installing, see docs/manual_steps.md\n"
            "    conda install -y -c salilab modeller\n"
            "    pip install csb-pandora\n"
            "Get a key from https://salilab.org/modeller/registration.html "
            f"({exc})"
        )
    except BaseException as exc:  # MODELLER may exit rather than raise on a bad key
        return False, f"MODELLER is installed but could not be imported: {exc!r}"

    # Environ() prints a banner to stdout; swallow it so a readiness check stays quiet.
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            modeller.Environ()
    except BaseException as exc:
        return False, (
            f"MODELLER {modeller.__version__} is installed but its license is not valid "
            f"({exc!r}). Put your key on the `license` line of\n"
            "    /opt/anaconda3/lib/modeller-10.8/modlib/modeller/config.py\n"
            "replacing the XXXX placeholder. Note config.py takes precedence over the "
            "KEY_MODELLER environment variable, so exporting that after installing has "
            "no effect -- see docs/manual_steps.md."
        )

    try:
        import PANDORA  # noqa: F401
    except ImportError as exc:
        return False, f"MODELLER is usable but PANDORA is not (pip install csb-pandora): {exc}"

    return True, f"PANDORA and MODELLER {modeller.__version__} available, license valid"


# --- the PANDORA route ------------------------------------------------------------------
#
# MHC class II binds a 9-residue core with four anchor pockets, P1/P4/P6/P9, and these are
# their offsets inside that core (0-based). PANDORA wants them as 1-based positions in the
# *full* peptide, so the core's offset within the extended sequence is added back.
#
# Why four and not the five pockets named in the literature (P1/P4/P6/P7/P9): PANDORA's class
# II alignment uses exactly four anchors -- Align.py pads peptides to align the 2nd and 3rd --
# and P1/P4/P6/P9 are the canonical class II set. P7 matters chemically for DQ2.5 (it also
# takes a deamidation-derived negative charge) but it is not a geometric anchor here.
CLASS_II_ANCHOR_OFFSETS = (0, 3, 5, 8)

#: PANDORA allele names for HLA-DQ2.5 = DQA1*05:01 + DQB1*02:01, in (alpha, beta) order.
#: The spelling is IMGT's as PANDORA stores it: asterisk kept, four digits, no colon --
#: `HLA-DQA1*0501`, not `HLA-DQA1*05:01` and not `HLA-DQA10501`. That is not guessable, and
#: the first guess here was wrong; `check_allele_names()` caught it against the built database.
#: Keep that check in the path -- PANDORA does not reject an unknown allele name, it scores
#: every template and returns the best, so a misspelling yields a model built on some other
#: allele's structure rather than an error.
DQ25_ALLELES = ["HLA-DQA1*0501", "HLA-DQB1*0201"]


def class_ii_anchors(extended_sequence: str, core: str) -> list[int]:
    """1-based positions of the P1/P4/P6/P9 anchors within `extended_sequence`.

    Derived from where the core sits, never hardcoded: the panel's extended sequences put the
    core at different offsets (E01's core starts at position 3 of `LQPFPQPELPY`, E04's at
    position 1 of `PQPEQPFPWQPEQ`), so a fixed anchor list would be right for one row and
    quietly wrong for the next.
    """
    if len(core) != 9:
        raise ValueError(f"class II core must be 9 residues, got {len(core)}: {core!r}")
    start = extended_sequence.find(core)
    if start < 0:
        raise ValueError(
            f"core {core!r} is not a substring of {extended_sequence!r}. For a deamidated row "
            "check the extended sequence carries the same Q->E substitutions as the core -- "
            "that mismatch is what made E09 wrong in Checkpoint 1."
        )
    return [start + 1 + offset for offset in CLASS_II_ANCHOR_OFFSETS]


def allele_vocabulary(database) -> set[str]:
    """Every class II allele name appearing on any template in `database`."""
    known: set[str] = set()
    for template in database.MHCII_data.values():
        allele_type = template.allele_type
        known.update(allele_type if isinstance(allele_type, list) else [allele_type])
    return known


def check_allele_names(alleles: list[str], database) -> None:
    """Raise unless every name in `alleles` is a real allele name PANDORA knows.

    PANDORA does not complain about an unrecognized allele name, so a typo does not fail --
    it silently changes which templates are eligible. Hence the check.

    Validate against the **full** database, not a restricted subset: "this name is misspelled"
    and "no template in this subset carries this name" are different problems, and only the
    first is an error. `templates_for_alleles()` answers the second.
    """
    known = allele_vocabulary(database)
    missing = [name for name in alleles if name not in known]
    if not missing:
        return

    hints = {}
    for name in missing:
        stem = name.split("-")[-1][:6]
        hints[name] = sorted(k for k in known if stem in k)[:6]
    raise ValueError(
        f"allele name(s) not present in PANDORA's class II template database: {missing}.\n"
        f"Closest names actually in the database: {hints}\n"
        f"({len(known)} distinct class II alleles across {len(database.MHCII_data)} templates.)"
    )


def templates_for_alleles(alleles: list[str], database) -> list[str]:
    """Template IDs PANDORA will consider for `alleles`, reproducing its own filter.

    PANDORA's class II template selection (`Modelling_functions.find_template`) is a hard
    boolean filter followed by a peptide-alignment score:

        if any(x in template.allele_type for x in target.allele_type)

    Two consequences worth knowing, both non-obvious:

    1. It is `any`, not `all`. Naming both DQ2.5 chains does not narrow the set beyond what
       either name alone admits, and does not give a template matching both extra weight.
    2. PANDORA's beta-chain annotation is incomplete. Of the 12 templates carrying
       `HLA-DQA1*0501`, only 1S9V also records `HLA-DQB1*0201` -- the rest (4OZF/G/H, 4OZI,
       5KSA/B, 6MFF/G, 6U3M/N/O) are genuinely DQ2.5 structures whose beta allele PANDORA
       simply did not parse from the IMGT remarks.

    Together those make a trap: asking for `HLA-DQB1*0201` **alone** is a correct, specific
    name for DQ2.5 that would quietly reduce the candidate set to the single structure 1S9V,
    so every panel peptide would be modeled on it. Pairing it with `HLA-DQA1*0501` is what
    keeps the real candidate set in play -- see `DQ25_ALLELES`.
    """
    chosen = []
    for pdb_id, template in database.MHCII_data.items():
        allele_type = template.allele_type
        allele_type = allele_type if isinstance(allele_type, list) else [allele_type]
        if any(name in allele_type for name in alleles):
            chosen.append(pdb_id)
    return sorted(chosen)


def load_database(path: Path | None = None):
    """Load the pickled class II template database built by `grains.pandora_db`."""
    from . import pandora_compat  # noqa: F401  (must precede PANDORA)
    from . import pandora_db
    from PANDORA.Database import Database

    ready, detail = pandora_db.database_ready()
    if not ready:
        raise FileNotFoundError(detail)
    return Database.load(str(path or pandora_db.DATABASE_PATH))


def pandora_model(
    peptide: str,
    *,
    peptide_id: str,
    output_dir: Path,
    database=None,
    alleles: list[str] | None = None,
    allele_reference=None,
    core: str | None = None,
    anchors: list[int] | None = None,
    n_loop_models: int = 20,
    loop_refinement: str = "slow",
    n_jobs: int | None = None,
) -> dict:
    """Model one peptide in the DQ2.5 groove with PANDORA. Returns a result dict.

    This is the only route whose output may be PRODIGY-scored -- see `assert_scorable`.

    `core` or `anchors` must be given. Passing `core` is preferred: the anchors are then
    derived from it by `class_ii_anchors`, which fails loudly if the core and the extended
    sequence disagree.

    MHC class II targets must supply both chain sequences; PANDORA has no reference-sequence
    fallback for class II (it raises "Reference MHC II sequences have not been implemented
    yet"). They come from `load_hla_chains`, i.e. from 1S9V's SEQRES records, so a model is
    built on the same constructs the Checkpoint 2 RMSD was measured against.
    """
    from . import pandora_compat  # noqa: F401  (must precede PANDORA)
    from PANDORA.PMHC import PMHC
    from PANDORA.Pandora import Pandora

    if anchors is None:
        if core is None:
            raise ValueError("pass either `core` (preferred) or explicit `anchors`")
        anchors = class_ii_anchors(peptide, core)

    alleles = list(alleles or DQ25_ALLELES)
    database = database if database is not None else load_database()
    # Spelling is checked against the full vocabulary (`allele_reference` when the modeling
    # database is a deliberately restricted subset, as in the Checkpoint 3 gate).
    check_allele_names(alleles, allele_reference if allele_reference is not None else database)

    eligible = templates_for_alleles(alleles, database)
    if not eligible:
        raise ValueError(
            f"no template in this database carries any of {alleles}, so PANDORA has nothing "
            f"to model on ({len(database.MHCII_data)} templates available). PANDORA's own "
            "error for this is 'could not find any putative template'."
        )

    alpha, beta = load_hla_chains()
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    target = PMHC.Target(
        id=peptide_id,
        allele_type=alleles,
        peptide=peptide,
        MHC_class="II",
        M_chain_seq=alpha,
        N_chain_seq=beta,
        anchors=anchors,
    )

    case = Pandora.Pandora(target, database, output_dir=str(output_dir))
    started = time.time()
    case.model(
        n_loop_models=n_loop_models,
        loop_refinement=loop_refinement,
        n_jobs=n_jobs,
        benchmark=False,   # keeps Model.renumber() (and its stubbed Biopython calls) off the path
        verbose=True,
    )
    elapsed = time.time() - started

    # PANDORA writes its models into a per-run subdirectory; collect them by mtime so the
    # ranking PANDORA itself produced is preserved rather than guessed from filenames.
    models = sorted(output_dir.rglob("*.BL*.pdb"), key=lambda p: p.stat().st_mtime)
    if not models:
        models = sorted(output_dir.rglob("*.pdb"), key=lambda p: p.stat().st_mtime)

    return {
        "route": ROUTE_PANDORA,
        "peptide_id": peptide_id,
        "peptide": peptide,
        "alleles": alleles,
        "anchors": anchors,
        "n_models": len(models),
        "templates_eligible": len(eligible),
        "models": [str(p) for p in models],
        "template_id": getattr(getattr(case, "template", None), "id", None),
        "seconds": round(elapsed, 1),
        "output_dir": str(output_dir),
    }


def main() -> None:
    write_validation_submission()
    ok, detail = pandora_available()
    print()
    print(f"PANDORA route: {'ready' if ok else 'NOT ready'}")
    if not ok:
        print(f"  {detail}")


if __name__ == "__main__":
    main()
