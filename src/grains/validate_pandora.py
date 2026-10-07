"""Checkpoint 3 go/no-go gate: is PANDORA accurate enough for the full panel?

Models the validation peptide (E01, the molecule deposited in 1S9V) with PANDORA and measures
it against the same crystal structure, through the same apparatus and the same interpretive
scale Checkpoint 2 built in `grains.superpose`.

Why this is three runs and not one
----------------------------------
**1S9V is itself in PANDORA's template database** -- with the exact peptide `LQPFPQPELPY`. Ask
PANDORA to model that peptide with the full database and it picks 1S9V as the template and
threads the target onto it, which is close to copying the answer. The resulting RMSD would be
beautiful and would mean nothing. Checkpoint 2 asked the same question of AlphaFold (was 1S9V a
template?) and had to settle for measuring the consequences; here the template set is ours to
control, so the question can be answered directly by withholding it.

| Condition | Withheld | What it measures |
|---|---|---|
| `full_db` | nothing | The circular ceiling. Reported precisely so it can be dismissed. |
| `no_exact` | 1S9V | Can PANDORA rebuild the complex without the complex? |
| `no_core` | every template whose peptide contains the E01 core | The honest estimate for the panel. |

`no_core` is the condition the gate decision rests on, and it is the one that resembles the
real task. For a gluten epitope PANDORA will usually have a close DQ2.5 template, but for a
**control peptide from a safe grain** it will not -- no such complex has ever been
crystallized. If accuracy collapses once near-identical templates are withheld, then the
toxic-vs-safe comparison would be measuring template availability rather than binding, which
would be a confound severe enough to invalidate the design.

Every run also records the chosen template and the longest residue run it shares with the
target peptide, because that is the covariate which explains the RMSD.

Usage
-----
    PYTHONPATH=src python -m grains.validate_pandora

Output
------
    data/processed/rmsd_pandora.csv      one row per (condition, model, crystal copy)
    data/processed/pandora_gate.json     the gate decision inputs
"""

from __future__ import annotations

import copy
import csv
import json
import shutil
import time
from pathlib import Path

from . import config, modeling, pandora_db
from .superpose import (
    E01_CORE,
    FIELDNAMES,
    S9V_COPY1,
    S9V_COPY2,
    ComplexSpec,
    compare,
)

#: PANDORA writes MHC class II models with chain M = alpha, N = beta, P = peptide.
# Re-exported so callers reading this module don't have to know where the IMGT
# annotation lives.
is_tcr_bound = pandora_db.is_tcr_bound
tcr_bound = pandora_db.tcr_bound

PANDORA_CHAINS = ("M", "N", "P")

#: Final loop-refined models. `.IL*` is the initial loop model and `.B9999*` the raw homology
#: model; neither is the output PANDORA ranks.
MODEL_GLOB = "*.BL*.pdb"

WORK_DIR = config.RAW / "models" / "pandora_gate"
RMSD_PATH = config.PROCESSED / "rmsd_pandora.csv"
GATE_PATH = config.PROCESSED / "pandora_gate.json"

# From data/processed/rmsd_baseline.csv -- the scale Checkpoint 2 measured before any model
# existed. These are what make an RMSD here mean something.
FLOOR_A = 0.24       # two copies in one 1S9V asymmetric unit: refinement noise
CEILING_A = 1.16     # 1S9V vs 4OZI: two independent depositions of the same complex


def longest_shared_run(a: str, b: str) -> str:
    """The longest contiguous residue run common to two peptides.

    The covariate that explains a homology model's accuracy: a template sharing 11 residues
    with the target is nearly the answer, one sharing 3 is a genuine prediction.
    """
    best = ""
    for start in range(len(a)):
        for end in range(start + len(best) + 1, len(a) + 1):
            if a[start:end] in b:
                best = a[start:end]
            else:
                break
    return best


def templates_containing(database, motif: str) -> set[str]:
    """IDs of every template whose peptide contains `motif`."""
    return {
        pdb_id for pdb_id, template in database.MHCII_data.items()
        if motif in (template.peptide or "")
    }


def restricted(database, withhold: set[str]):
    """A copy of the database with `withhold` removed from the class II template set."""
    subset = copy.copy(database)
    subset.MHCII_data = {
        pdb_id: template for pdb_id, template in database.MHCII_data.items()
        if pdb_id not in withhold
    }
    return subset


def conditions(database) -> dict[str, dict]:
    """The three template-withholding conditions, built against the real database."""
    core_bearing = templates_containing(database, E01_CORE)
    return {
        "full_db": {
            "withhold": set(),
            "rationale": "nothing withheld; 1S9V is a template, so this is circular",
            "circular": True,
        },
        "no_exact": {
            "withhold": {"1S9V"},
            "rationale": "the exact deposited complex withheld",
            "circular": False,
        },
        "no_core": {
            "withhold": set(core_bearing),
            "rationale": (f"every template whose peptide contains the E01 core {E01_CORE} "
                          f"withheld ({len(core_bearing)}: "
                          f"{', '.join(sorted(core_bearing))})"),
            "circular": False,
        },
        # The controlled comparison for the conformational-state effect. `no_exact` fell back
        # to 4OZI, which shares 10 of 11 residues but is TCR-bound; this condition withholds
        # 1S9V *and* every ternary template, so the best remaining candidate has comparable
        # sequence similarity in the binary state. Similarity held roughly constant, state
        # varied -- which is what isolates the effect rather than inferring it.
        "no_exact_binary": {
            "withhold": {"1S9V"} | pandora_db.tcr_bound(database),
            "rationale": ("1S9V and every TCR-bound (ternary) template withheld; isolates "
                          "conformational state from sequence similarity"),
            "circular": False,
        },
    }


def run_condition(name: str, spec: dict, database, *, n_loop_models: int) -> dict:
    """Model E01 under one condition and measure it against both 1S9V crystal copies."""
    print(f"\n{'=' * 78}\n{name}: {spec['rationale']}\n{'=' * 78}", flush=True)

    out_dir = WORK_DIR / name
    if out_dir.exists():
        shutil.rmtree(out_dir)          # a stale run would be silently picked up by the glob

    subset = restricted(database, spec["withhold"])
    started = time.time()
    result = modeling.pandora_model(
        modeling.VALIDATION_PEPTIDE,
        peptide_id=f"E01_{name}",
        core=E01_CORE,
        output_dir=out_dir,
        database=subset,
        allele_reference=database,   # spelling checked against the full vocabulary
        n_loop_models=n_loop_models,
    )
    elapsed = time.time() - started

    template_id = result.get("template_id")
    template = subset.MHCII_data.get(template_id) if template_id else None
    template_peptide = (template.peptide or "") if template else ""
    shared = longest_shared_run(modeling.VALIDATION_PEPTIDE, template_peptide)

    model_paths = sorted(out_dir.rglob(MODEL_GLOB))
    rows: list[dict] = []
    for path in model_paths:
        model = ComplexSpec(f"{name}:{path.stem.split('.')[-1]}", path, *PANDORA_CHAINS)
        for crystal in (S9V_COPY1, S9V_COPY2):
            row = compare(model, crystal)
            row["note"] = f"{name}; template {template_id}; shared run {len(shared)} aa"
            rows.append(row)

    cores = [row["core_rmsd_A"] for row in rows]
    summary = {
        "condition": name,
        "rationale": spec["rationale"],
        "circular": spec["circular"],
        "withheld": sorted(spec["withhold"]),
        "templates_available": len(subset.MHCII_data),
        "template_chosen": template_id,
        "template_peptide": template_peptide,
        "shared_run": shared,
        "shared_run_len": len(shared),
        "n_models": len(model_paths),
        "core_rmsd_min_A": round(min(cores), 2) if cores else None,
        "core_rmsd_max_A": round(max(cores), 2) if cores else None,
        "seconds": round(elapsed, 1),
    }

    print(f"\n  template chosen   {template_id}  ({template_peptide})")
    print(f"  shared run        {shared!r} ({len(shared)} aa of "
          f"{len(modeling.VALIDATION_PEPTIDE)})")
    if cores:
        print(f"  core RMSD         {min(cores):.2f}-{max(cores):.2f} Å "
              f"(floor {FLOOR_A}, ceiling {CEILING_A})")
    print(f"  wall clock        {elapsed:.0f} s for {len(model_paths)} models")
    return {"summary": summary, "rows": rows}


def run(n_loop_models: int = 20) -> dict:
    config.ensure_dirs()
    WORK_DIR.mkdir(parents=True, exist_ok=True)

    database = modeling.load_database()
    modeling.check_allele_names(modeling.DQ25_ALLELES, database)
    print(f"template database: {len(database.MHCII_data)} MHC class II structures")

    specs = conditions(database)
    results = [run_condition(name, spec, database, n_loop_models=n_loop_models)
               for name, spec in specs.items()]

    rows = [row for result in results for row in result["rows"]]
    with RMSD_PATH.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDNAMES)  # already ends with "note"
        writer.writeheader()
        writer.writerows(rows)

    summaries = [result["summary"] for result in results]
    honest = next(s for s in summaries if s["condition"] == "no_core")
    report = {
        "run_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "peptide": modeling.VALIDATION_PEPTIDE,
        "core": E01_CORE,
        "n_loop_models": n_loop_models,
        "reference_floor_A": FLOOR_A,
        "reference_ceiling_A": CEILING_A,
        "conditions": summaries,
        "gate_metric": "no_core core_rmsd_max_A",
        "gate_value_A": honest["core_rmsd_max_A"],
        "gate_pass": (honest["core_rmsd_max_A"] is not None
                      and honest["core_rmsd_max_A"] <= CEILING_A),
    }
    GATE_PATH.write_text(json.dumps(report, indent=2) + "\n")

    print(f"\n{'=' * 78}\nsummary\n{'=' * 78}")
    print(f"{'condition':<17}{'template':<10}{'shared':>7}{'core RMSD (Å)':>16}{'sec':>8}")
    for summary in summaries:
        span = (f"{summary['core_rmsd_min_A']:.2f}-{summary['core_rmsd_max_A']:.2f}"
                if summary["core_rmsd_min_A"] is not None else "n/a")
        flag = "  (circular)" if summary["circular"] else ""
        print(f"{summary['condition']:<17}{str(summary['template_chosen']):<10}"
              f"{summary['shared_run_len']:>7}{span:>16}{summary['seconds']:>8.0f}{flag}")
    print(f"\nreference: floor {FLOOR_A} Å (one crystal, two copies), "
          f"ceiling {CEILING_A} Å (independent crystals)")
    print(f"gate metric: no_core worst core RMSD = {report['gate_value_A']} Å -> "
          f"{'PASS' if report['gate_pass'] else 'FAIL'}")
    print(f"\nwrote {RMSD_PATH}\nwrote {GATE_PATH}")
    return report


if __name__ == "__main__":
    run()
