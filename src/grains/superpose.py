"""Superpose two peptide-HLA complexes and report backbone RMSD.

Checkpoint 2's question is "can the modeling approach reproduce a known complex?", which
is a number: the RMSD between a model and the 1S9V crystal structure. But an RMSD on its
own is uninterpretable -- is 1.5 A good? Good compared to what? A model can only be asked
to agree with a crystal structure about as well as crystal structures of the same thing
agree with each other, so this module first measures that experimental noise floor from
structures alone, before any model exists:

  * 1S9V copy 1 vs copy 2   -- two copies in one asymmetric unit, same crystal, same
                               refinement. The tightest floor available: differences here
                               are refinement noise, nothing else.
  * 4OZI copy 1 vs copy 2   -- the same, in the lower-resolution (3.2 A) ternary crystal.
  * 1S9V      vs 4OZI       -- independent depositions, different crystal forms, different
                               labs, and binary (no TCR) vs TCR-bound. Both carry the E01
                               (DQ2.5-glia-a1a) peptide, so this is the honest ceiling on
                               agreement: no model should be expected to beat it.

A model landing inside that spread has reproduced the complex. One landing well outside it
has not, however small the number looks in isolation.

Protocol (the part that is easy to get wrong)
---------------------------------------------
Fit on the **MHC groove** (alpha + beta chain CA atoms), then measure the peptide RMSD
*where it lands* under that fit. Do not re-fit on the peptide: a peptide-to-peptide fit
reports how similar the peptide's internal shape is while discarding whether it sits in
the same place in the groove, which is the entire question. The MHC fit is the common
reference frame; the peptide is the measurement.

Residue correspondence comes from **sequence alignment**, never residue numbers. Numbering
is not comparable across these depositions: the alpha chains carry insertion codes, and
the peptide chains are numbered 1..11 in 1S9V but 2..12 in 4OZI for different sequences
(4OZI's peptide also carries a C-terminal PGS construct tag that is not native gliadin).
Aligning sequences and keeping only identically-matched columns is robust to all of that.

Usage
-----
    PYTHONPATH=src python -m grains.superpose

Output
------
    data/processed/rmsd_baseline.csv   one row per comparison
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from Bio.Align import PairwiseAligner
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import is_aa
from Bio.PDB.Superimposer import Superimposer

from . import config

PDB_DIR = config.RAW / "pdb"

# The E01 (DQ2.5-glia-a1a) 9mer core, present in both validation structures' peptides.
# Peptide RMSD is reported over this window as well as over the whole resolved peptide,
# since the core is what occupies the anchor pockets and what the scoring step cares about.
E01_CORE = "PFPQPELPY"

# Three-letter -> one-letter, for reading sequences off the resolved ATOM records.
THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}


@dataclass(frozen=True)
class ComplexSpec:
    """One peptide-MHC complex: which file, and which chains play which role."""

    label: str
    pdb_path: Path
    alpha_chain: str
    beta_chain: str
    peptide_chain: str


def resolved_residues(structure, chain_id: str) -> tuple[list, str]:
    """Residues of `chain_id` that have a CA atom, in order, plus their sequence.

    Only standard amino acids are kept, so waters, ions and any hetero-residues are
    dropped. Returns the residue objects alongside the one-letter sequence so an index
    into the sequence maps straight back to a residue.
    """
    model = next(structure.get_models())
    residues = [
        residue
        for residue in model[chain_id]
        if is_aa(residue, standard=True) and residue.has_id("CA")
    ]
    sequence = "".join(THREE_TO_ONE[r.get_resname()] for r in residues)
    return residues, sequence


def matched_pairs(seq_a: str, seq_b: str) -> list[tuple[int, int]]:
    """Index pairs (i in seq_a, j in seq_b) for identically-matched aligned residues.

    A global alignment, then only the columns where both sequences have a residue and
    those residues are the same. Mismatched columns are discarded rather than paired:
    these are near-identical constructs, so a mismatch means a construct difference or a
    misalignment, and either way such a column should not be driving a structural fit.
    """
    aligner = PairwiseAligner(
        mode="global",
        match_score=1,
        mismatch_score=-1,
        open_gap_score=-5,
        extend_gap_score=-0.5,
    )
    alignment = aligner.align(seq_a, seq_b)[0]
    pairs: list[tuple[int, int]] = []
    for (a_start, a_end), (b_start, b_end) in zip(*alignment.aligned):
        for offset in range(a_end - a_start):
            i, j = a_start + offset, b_start + offset
            if seq_a[i] == seq_b[j]:
                pairs.append((i, j))
    return pairs


def ca_coords(residues: list, indices: list[int]) -> np.ndarray:
    return np.array([residues[i]["CA"].get_coord() for i in indices], dtype=float)


# A peptide bound in the groove has its closest CA within ~5 A of both MHC chains (measured:
# 4.8-5.4 A across both validation structures). A peptide belonging to a different copy in
# the asymmetric unit sits 16-30 A away. 15 A separates those cleanly with room to spare.
MAX_BOUND_CA_DISTANCE_A = 15.0


def verify_pairing(structure, spec: "ComplexSpec") -> dict[str, float]:
    """Check the spec's peptide chain really is bound to its alpha/beta chains.

    Chain-ID order does NOT imply which chains form a complex when the asymmetric unit
    holds more than one copy. In 4OZI it is actively misleading: peptide chain J binds
    MHC A/B and peptide chain I binds C/D, the opposite of the lettering (1S9V, by
    contrast, pairs C with A/B and F with D/E as the order suggests). Pairing A/B with I
    on the strength of the lettering superposes the groove correctly and then measures a
    peptide from the neighbouring complex ~62 A away -- a wrong answer that looks like a
    computed result rather than an error, which is why this is a hard failure and not a
    warning.
    """
    peptide_res, _ = resolved_residues(structure, spec.peptide_chain)
    peptide_xyz = np.array([r["CA"].get_coord() for r in peptide_res], dtype=float)

    distances: dict[str, float] = {}
    for role, chain_id in (("alpha", spec.alpha_chain), ("beta", spec.beta_chain)):
        mhc_res, _ = resolved_residues(structure, chain_id)
        mhc_xyz = np.array([r["CA"].get_coord() for r in mhc_res], dtype=float)
        pairwise = np.linalg.norm(peptide_xyz[:, None, :] - mhc_xyz[None, :, :], axis=-1)
        distances[role] = float(pairwise.min())

    too_far = {role: d for role, d in distances.items() if d > MAX_BOUND_CA_DISTANCE_A}
    if too_far:
        detail = ", ".join(f"{role} chain min CA distance {d:.1f} A" for role, d in distances.items())
        raise ValueError(
            f"{spec.label}: peptide chain {spec.peptide_chain} is not bound to "
            f"alpha={spec.alpha_chain}/beta={spec.beta_chain} ({detail}; "
            f"threshold {MAX_BOUND_CA_DISTANCE_A} A). The peptide almost certainly belongs "
            "to a different copy in the asymmetric unit -- check which chains actually "
            "contact each other instead of trusting chain-ID order."
        )
    return distances


def rmsd(coords_a: np.ndarray, coords_b: np.ndarray) -> float:
    """Plain RMSD between two already-superposed, already-paired coordinate sets."""
    diff = coords_a - coords_b
    return float(np.sqrt((diff * diff).sum() / len(diff)))


def core_window(sequence: str, core: str) -> list[int] | None:
    """Indices of `core` within `sequence`, or None when it isn't present."""
    start = sequence.find(core)
    if start < 0:
        return None
    return list(range(start, start + len(core)))


def parser_for(path: Path):
    """Biopython parser matching the file extension.

    AlphaFold Server returns mmCIF; the crystal structures here are legacy PDB. Both parse
    into the same Structure object, so everything downstream is format-agnostic.
    """
    if path.suffix.lower() in (".cif", ".mmcif"):
        from Bio.PDB import MMCIFParser

        return MMCIFParser(QUIET=True)
    return PDBParser(QUIET=True)


def compare(mobile: ComplexSpec, target: ComplexSpec, core: str = E01_CORE) -> dict:
    """Fit `mobile` onto `target` by their MHC chains; report where the peptide lands.

    Returns the fit RMSD over the groove, plus peptide RMSD over the whole matched
    peptide and over the `core` window, all in Angstroms, with the residue counts each
    number is based on (an RMSD without its N is not reportable).
    """
    mobile_structure = parser_for(mobile.pdb_path).get_structure(
        mobile.label, mobile.pdb_path)
    target_structure = parser_for(target.pdb_path).get_structure(
        target.label, target.pdb_path)

    # Refuse to report an RMSD for chains that aren't actually a complex (see verify_pairing).
    verify_pairing(mobile_structure, mobile)
    verify_pairing(target_structure, target)

    # --- fit on the groove: alpha and beta chains together -------------------------
    fixed_atoms, moving_atoms = [], []
    for mobile_chain, target_chain in (
        (mobile.alpha_chain, target.alpha_chain),
        (mobile.beta_chain, target.beta_chain),
    ):
        mobile_res, mobile_seq = resolved_residues(mobile_structure, mobile_chain)
        target_res, target_seq = resolved_residues(target_structure, target_chain)
        for i, j in matched_pairs(mobile_seq, target_seq):
            moving_atoms.append(mobile_res[i]["CA"])
            fixed_atoms.append(target_res[j]["CA"])

    superimposer = Superimposer()
    superimposer.set_atoms(fixed_atoms, moving_atoms)
    rotation, translation = superimposer.rotran

    # --- measure the peptide under that fit ----------------------------------------
    # The transform is applied to copied coordinates rather than with Superimposer.apply(),
    # which mutates the parsed structure in place; keeping the structure pristine means a
    # caller can run several comparisons off one parse without silent carry-over.
    mobile_pep, mobile_pep_seq = resolved_residues(mobile_structure, mobile.peptide_chain)
    target_pep, target_pep_seq = resolved_residues(target_structure, target.peptide_chain)
    pep_pairs = matched_pairs(mobile_pep_seq, target_pep_seq)

    mobile_idx = [i for i, _ in pep_pairs]
    target_idx = [j for _, j in pep_pairs]
    mobile_pep_xyz = ca_coords(mobile_pep, mobile_idx) @ rotation + translation
    target_pep_xyz = ca_coords(target_pep, target_idx)
    peptide_rmsd = rmsd(mobile_pep_xyz, target_pep_xyz)

    # Core RMSD: restrict to pairs where BOTH sides fall inside their own core window,
    # so the number is over the core in both structures, not one structure's core
    # compared against whatever the other happened to align there.
    mobile_core = core_window(mobile_pep_seq, core)
    target_core = core_window(target_pep_seq, core)
    if mobile_core is None or target_core is None:
        core_rmsd: float | None = None
        n_core = 0
    else:
        core_pairs = [
            (i, j) for i, j in pep_pairs if i in set(mobile_core) and j in set(target_core)
        ]
        if core_pairs:
            core_mobile = ca_coords(mobile_pep, [i for i, _ in core_pairs]) @ rotation + translation
            core_target = ca_coords(target_pep, [j for _, j in core_pairs])
            core_rmsd = rmsd(core_mobile, core_target)
            n_core = len(core_pairs)
        else:
            core_rmsd, n_core = None, 0

    return {
        "comparison": f"{mobile.label} -> {target.label}",
        "mobile": mobile.label,
        "target": target.label,
        "n_mhc_ca": len(fixed_atoms),
        "mhc_fit_rmsd_A": round(superimposer.rms, 3),
        "n_peptide_ca": len(pep_pairs),
        "peptide_rmsd_A": round(peptide_rmsd, 3),
        "n_core_ca": n_core,
        "core_rmsd_A": None if core_rmsd is None else round(core_rmsd, 3),
        "mobile_peptide_seq": mobile_pep_seq,
        "target_peptide_seq": target_pep_seq,
    }


# --- the crystal-vs-crystal baseline -----------------------------------------------
# 1S9V: binary DQ2.5 + E01, 2.2 A, two copies in the ASU. Peptide C binds A/B, F binds D/E.
# 4OZI: ternary DQ2.5 + E01 + TCR s2, 3.2 A, two copies. Peptide J binds A/B and I binds C/D
#       -- the reverse of what the chain lettering suggests; confirmed by contact distance,
#       and enforced by verify_pairing() so the swap can't be reintroduced unnoticed.
S9V_COPY1 = ComplexSpec("1S9V_ABC", PDB_DIR / "1S9V.pdb", "A", "B", "C")
S9V_COPY2 = ComplexSpec("1S9V_DEF", PDB_DIR / "1S9V.pdb", "D", "E", "F")
OZI_COPY1 = ComplexSpec("4OZI_ABJ", PDB_DIR / "4OZI.pdb", "A", "B", "J")
OZI_COPY2 = ComplexSpec("4OZI_CDI", PDB_DIR / "4OZI.pdb", "C", "D", "I")

BASELINE_COMPARISONS = [
    (S9V_COPY1, S9V_COPY2, "same crystal (1S9V ASU copies) -- refinement noise floor"),
    (OZI_COPY1, OZI_COPY2, "same crystal (4OZI ASU copies) at 3.2 A"),
    (S9V_COPY1, OZI_COPY1, "independent depositions; binary vs TCR-bound"),
    (S9V_COPY2, OZI_COPY2, "independent depositions, second copy of each"),
]

FIELDNAMES = [
    "comparison", "mobile", "target", "n_mhc_ca", "mhc_fit_rmsd_A",
    "n_peptide_ca", "peptide_rmsd_A", "n_core_ca", "core_rmsd_A",
    "mobile_peptide_seq", "target_peptide_seq", "note",
]


AF_MODEL_DIR = config.RAW / "models" / "fold_dq2_5_e01_validation_1s9v"
AF_MODEL_GLOB = "fold_dq2_5_e01_validation_1s9v_model_*.cif"

# Structures AlphaFold Server used as templates for this job, read from the returned
# templates/ directory. 1S9V is deliberately absent, which is what keeps the validation from
# being circular; 4OZI is present for the peptide chain, which is a caveat to state rather
# than hide -- see validate_models().
AF_TEMPLATES = {"6px6", "4ozi", "4d8p", "6dig", "6u3m", "3pl6", "1uvq", "1iak", "4grl"}


def validate_models(model_dir: Path | None = None) -> list[dict]:
    """Compare every AlphaFold model against both crystal copies of 1S9V and 4OZI.

    Checkpoint 2's deliverable. Each model is compared with the validation target (1S9V,
    *not* an AlphaFold template) and with 4OZI (which *was* a template for the peptide
    chain). Reporting both is the point: if the model had simply recalled its template it
    would sit at 4OZI's own distance from 1S9V (~1.1 A core) and be much closer to 4OZI than
    to 1S9V. Being roughly equidistant from the two, and closer to each than they are to each
    other, is the signature of a consensus prediction instead.
    """
    model_dir = model_dir or AF_MODEL_DIR
    paths = sorted(model_dir.glob(AF_MODEL_GLOB))
    if not paths:
        raise FileNotFoundError(
            f"no AlphaFold models in {model_dir} (expected {AF_MODEL_GLOB}) -- "
            "see docs/manual_steps.md for the submission and download steps"
        )

    targets = [
        (S9V_COPY1, "validation target (not an AlphaFold template)"),
        (S9V_COPY2, "validation target (not an AlphaFold template)"),
        (OZI_COPY1, "WAS an AlphaFold template for the peptide chain"),
        (OZI_COPY2, "WAS an AlphaFold template for the peptide chain"),
    ]

    rows: list[dict] = []
    for path in paths:
        label = path.stem.replace("fold_dq2_5_e01_validation_1s9v_", "AF_")
        model = ComplexSpec(label, path, "A", "B", "C")
        for target, note in targets:
            row = compare(model, target)
            row["note"] = note
            rows.append(row)

    out_path = config.PROCESSED / "rmsd_validation.csv"
    with out_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDNAMES)
        writer.writeheader()
        writer.writerows(rows)

    print(f"{'model':<14}{'target':<12}{'groove':>8}{'peptide':>9}{'core':>8}")
    for row in rows:
        print(f"{row['mobile']:<14}{row['target']:<12}"
              f"{row['mhc_fit_rmsd_A']:>8.2f}{row['peptide_rmsd_A']:>9.2f}"
              f"{row['core_rmsd_A']:>8.2f}")

    against_1s9v = [r["core_rmsd_A"] for r in rows if r["target"].startswith("1S9V")]
    against_4ozi = [r["core_rmsd_A"] for r in rows if r["target"].startswith("4OZI")]
    print(f"\nwrote {out_path}")
    print(f"core RMSD vs 1S9V (unseen): {min(against_1s9v):.2f}-{max(against_1s9v):.2f} A")
    print(f"core RMSD vs 4OZI (template): {min(against_4ozi):.2f}-{max(against_4ozi):.2f} A")
    return rows


def run() -> list[dict]:
    config.ensure_dirs()
    missing = [
        spec.pdb_path
        for spec in (S9V_COPY1, OZI_COPY1)
        if not spec.pdb_path.exists()
    ]
    if missing:
        raise FileNotFoundError(
            "missing PDB file(s): "
            + ", ".join(str(p) for p in missing)
            + " -- download with e.g. "
            "curl https://files.rcsb.org/download/4OZI.pdb -o data/raw/pdb/4OZI.pdb"
        )

    rows = []
    for mobile, target, note in BASELINE_COMPARISONS:
        row = compare(mobile, target)
        row["note"] = note
        rows.append(row)
        print(f"{row['comparison']:>22}  "
              f"groove {row['mhc_fit_rmsd_A']:>5} A (n={row['n_mhc_ca']:>3})  "
              f"peptide {row['peptide_rmsd_A']:>5} A (n={row['n_peptide_ca']:>2})  "
              f"core {row['core_rmsd_A']:>5} A (n={row['n_core_ca']})")

    out_path = config.PROCESSED / "rmsd_baseline.csv"
    with out_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDNAMES)
        writer.writeheader()
        writer.writerows(rows)
    print(f"\nwrote {out_path}")

    peptide_values = [r["peptide_rmsd_A"] for r in rows]
    core_values = [r["core_rmsd_A"] for r in rows if r["core_rmsd_A"] is not None]
    print(f"peptide CA RMSD across comparisons: "
          f"{min(peptide_values):.2f}-{max(peptide_values):.2f} A")
    if core_values:
        print(f"E01 core CA RMSD across comparisons: "
              f"{min(core_values):.2f}-{max(core_values):.2f} A")
    return rows


if __name__ == "__main__":
    run()
