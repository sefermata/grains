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


def main() -> None:
    write_validation_submission()
    ok, detail = pandora_available()
    print()
    print(f"PANDORA route: {'ready' if ok else 'NOT ready'}")
    if not ok:
        print(f"  {detail}")


if __name__ == "__main__":
    main()
