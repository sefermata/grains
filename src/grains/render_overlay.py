"""Checkpoint 2 figure: AlphaFold model vs the 1S9V crystal structure.

Two panels, because the 3D overlay alone can't answer the question it raises.

Left -- the peptide backbones superposed in the groove, after fitting on the MHC chains.
Shows *that* the model puts the peptide in the right place.

Right -- per-residue CA deviation, which shows *where* it doesn't, against the two
reference bands measured in `grains.superpose`: the within-crystal refinement floor (two
copies of 1S9V, 0.24 A on the core) and the cross-deposition ceiling (1S9V vs 4OZI, ~1.16 A).
A deviation inside the lower band is as close as two refinements of the same crystal; inside
the upper band is as close as two independent crystals of the same complex. That's what makes
the number interpretable rather than just small.

The per-residue view also carries a real structural expectation: for MHC class II the 9-mer
core sits in the groove while the flanking residues extend out of it into solvent, so the
flanks are less constrained and should deviate more. Seeing that pattern is a sanity check on
the model; seeing the reverse would be a warning.

Colors are categorical slots 1 and 2 of the validated default palette, in fixed order
(crystal = blue, model = orange). CVD separation was computed, not eyeballed: normal-vision
dE 33.6/31.8 and protan/deutan/tritan dE 24.6-43.7, all above threshold in both modes. Light
and dark are rendered as separate files with their own steps -- a PNG can't respond to
prefers-color-scheme, so the dark variant is selected rather than an inverted light one.

Usage
-----
    PYTHONPATH=src python -m grains.render_overlay

Output
------
    results/figures/af_vs_1S9V_overlay.png
    results/figures/af_vs_1S9V_overlay_dark.png
"""

from __future__ import annotations

import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from . import config  # noqa: E402
from .superpose import (  # noqa: E402
    AF_MODEL_DIR,
    E01_CORE,
    S9V_COPY1,
    ComplexSpec,
    ca_coords,
    matched_pairs,
    parser_for,
    resolved_residues,
    verify_pairing,
)
from Bio.PDB.Superimposer import Superimposer  # noqa: E402

# Reference bands from data/processed/rmsd_baseline.csv (E01 core CA RMSD).
FLOOR_A = 0.24       # two copies in one 1S9V asymmetric unit -- refinement noise
CEILING_A = 1.16     # 1S9V vs 4OZI -- independent depositions of the same complex

THEMES = {
    "light": dict(
        surface="#fcfcfb", ink="#0b0b0b", ink2="#52514e", muted="#8a8a84",
        crystal="#2a78d6", model="#eb6834", band="#d8d8d2", suffix="",
    ),
    "dark": dict(
        surface="#1a1a19", ink="#ffffff", ink2="#c3c2b7", muted="#7e7d76",
        crystal="#3987e5", model="#d95926", band="#34332f", suffix="_dark",
    ),
}


def fit_model_to_crystal(model: ComplexSpec, crystal: ComplexSpec):
    """Fit on the groove; return the model's peptide CAs in the crystal's frame.

    Same protocol as `grains.superpose.compare` -- fit on MHC, measure the peptide where it
    lands, never re-fit on the peptide.
    """
    model_structure = parser_for(model.pdb_path).get_structure(model.label, model.pdb_path)
    crystal_structure = parser_for(crystal.pdb_path).get_structure(
        crystal.label, crystal.pdb_path)
    verify_pairing(model_structure, model)
    verify_pairing(crystal_structure, crystal)

    fixed, moving = [], []
    for model_chain, crystal_chain in (
        (model.alpha_chain, crystal.alpha_chain),
        (model.beta_chain, crystal.beta_chain),
    ):
        model_res, model_seq = resolved_residues(model_structure, model_chain)
        crystal_res, crystal_seq = resolved_residues(crystal_structure, crystal_chain)
        for i, j in matched_pairs(model_seq, crystal_seq):
            moving.append(model_res[i]["CA"])
            fixed.append(crystal_res[j]["CA"])

    superimposer = Superimposer()
    superimposer.set_atoms(fixed, moving)
    rotation, translation = superimposer.rotran

    model_pep, model_pep_seq = resolved_residues(model_structure, model.peptide_chain)
    crystal_pep, crystal_pep_seq = resolved_residues(
        crystal_structure, crystal.peptide_chain)
    pairs = matched_pairs(model_pep_seq, crystal_pep_seq)

    model_xyz = ca_coords(model_pep, [i for i, _ in pairs]) @ rotation + translation
    crystal_xyz = ca_coords(crystal_pep, [j for _, j in pairs])
    labels = [crystal_pep_seq[j] for _, j in pairs]

    # Groove chains, in the crystal's frame, as recessive context for the 3D panel.
    groove = []
    for chain_id in (crystal.alpha_chain, crystal.beta_chain):
        residues, _ = resolved_residues(crystal_structure, chain_id)
        groove.append(np.array([r["CA"].get_coord() for r in residues], dtype=float))

    return model_xyz, crystal_xyz, labels, groove, superimposer.rms, crystal_pep_seq


def render(theme_name: str, model_path: Path) -> Path:
    theme = THEMES[theme_name]
    model = ComplexSpec("AlphaFold model", model_path, "A", "B", "C")
    model_xyz, crystal_xyz, labels, groove, groove_rms, pep_seq = fit_model_to_crystal(
        model, S9V_COPY1)

    deviations = np.linalg.norm(model_xyz - crystal_xyz, axis=1)
    core_start = pep_seq.find(E01_CORE)
    core_range = range(core_start, core_start + len(E01_CORE)) if core_start >= 0 else range(0)
    core_rmsd = float(np.sqrt((deviations[list(core_range)] ** 2).mean())) if core_start >= 0 else float("nan")

    fig = plt.figure(figsize=(12.0, 5.2), facecolor=theme["surface"])
    grid = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.25], wspace=0.14)

    # --- left: 3D overlay ---------------------------------------------------------
    ax3d = fig.add_subplot(grid[0], projection="3d", facecolor=theme["surface"])
    # Show only groove residues near the peptide: the full CA trace of both chains is a
    # tangle at this scale and hides the thing the panel exists to show.
    centroid = crystal_xyz.mean(axis=0)
    for trace in groove:
        near = trace[np.linalg.norm(trace - centroid, axis=1) < 16.0]
        if len(near):
            ax3d.scatter(*near.T, color=theme["muted"], s=5.0, alpha=0.42,
                         depthshade=False, linewidths=0)
    ax3d.plot(*crystal_xyz.T, color=theme["crystal"], lw=2.4,
              marker="o", markersize=4.5, label="1S9V crystal")
    ax3d.plot(*model_xyz.T, color=theme["model"], lw=2.4,
              marker="o", markersize=4.5, label="AlphaFold model")
    ax3d.set_title("Peptide in the DQ2.5 groove\nfit on MHC chains, peptide measured",
                   color=theme["ink"], fontsize=10, pad=2)
    ax3d.set_box_aspect((1, 1, 1), zoom=1.25)
    # No axes on this panel. It answers "do the two backbones coincide?", which is
    # qualitative; absolute Angstrom coordinates are noise here and their tick labels
    # collide with the quantitative panel beside it. The numbers live on the right.
    ax3d.set_axis_off()
    ax3d.view_init(elev=22, azim=-62)
    legend = ax3d.legend(loc="upper left", fontsize=8, frameon=False)
    for text in legend.get_texts():
        text.set_color(theme["ink2"])

    # --- right: per-residue deviation ---------------------------------------------
    ax = fig.add_subplot(grid[1], facecolor=theme["surface"])
    positions = np.arange(len(deviations))
    top = max(float(deviations.max()), CEILING_A) * 1.30
    ax.set_xlim(-0.8, len(deviations) - 0.2)
    ax.set_ylim(0, top)

    # The core sits in the groove; the flanks hang out of it. One light band, drawn
    # first and behind everything, marks that boundary.
    if core_start >= 0:
        ax.axvspan(core_start - 0.5, core_start + len(E01_CORE) - 0.5,
                   color=theme["crystal"], alpha=0.07, lw=0, zorder=0)
        ax.text(core_start + len(E01_CORE) / 2 - 0.5, top * 0.985,
                "9-mer core — sits in the groove", ha="center", va="top",
                fontsize=8, color=theme["ink2"], zorder=4)

    # Reference levels as lines, not nested translucent fills: filled bands stack into
    # unreadable mud where they cross the core shading, and the thing being compared is
    # a level, which a line states directly.
    for level, text, side in ((CEILING_A, f"independent crystals  {CEILING_A} Å", "bottom"),
                              (FLOOR_A, f"one crystal, two copies  {FLOOR_A} Å", "top")):
        ax.axhline(level, color=theme["muted"], lw=1.0, ls=(0, (5, 3)), zorder=1)
        ax.text(-0.6, level, f" {text}", ha="left", va=side,
                fontsize=7.5, color=theme["muted"], zorder=4)

    ax.plot(positions, deviations, color=theme["model"], lw=2.0,
            marker="o", markersize=7, markeredgecolor=theme["surface"],
            markeredgewidth=2, zorder=3, label="AlphaFold model vs 1S9V")

    ax.set_xticks(positions)
    ax.set_xticklabels(list(labels), fontsize=9, color=theme["ink2"])
    ax.set_xlabel("peptide residue", color=theme["ink2"], fontsize=9)
    ax.set_ylabel("CA deviation from crystal (Å)", color=theme["ink2"], fontsize=9)
    ax.set_title(f"Per-residue deviation — core RMSD {core_rmsd:.2f} Å",
                 color=theme["ink"], fontsize=10, pad=8)
    ax.tick_params(colors=theme["muted"], labelsize=8)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(theme["muted"])
    ax.grid(axis="y", color=theme["muted"], alpha=0.15, lw=0.6)
    ax.set_axisbelow(True)

    fig.suptitle(
        "AlphaFold Server model reproduces the DQ2.5 + glia-α1a crystal structure",
        color=theme["ink"], fontsize=12, y=0.99)
    fig.text(0.5, 0.015,
             f"groove fit {groove_rms:.2f} Å over {len(groove[0]) + len(groove[1])} CA  ·  "
             "1S9V was not an AlphaFold template  ·  geometric comparison only, not a binding score",
             ha="center", fontsize=7.5, color=theme["muted"])
    fig.tight_layout(rect=(0, 0.05, 1, 0.94))

    out_dir = config.RESULTS / "figures"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"af_vs_1S9V_overlay{theme['suffix']}.png"
    fig.savefig(out_path, dpi=200, facecolor=theme["surface"])
    plt.close(fig)
    print(f"wrote {out_path}  (core RMSD {core_rmsd:.2f} Å, groove fit {groove_rms:.2f} Å)")
    return out_path


def run() -> list[Path]:
    warnings.filterwarnings("ignore")
    models = sorted(AF_MODEL_DIR.glob("*_model_0.cif"))
    if not models:
        raise FileNotFoundError(
            f"no AlphaFold model_0 in {AF_MODEL_DIR} -- see docs/manual_steps.md")
    return [render(name, models[0]) for name in THEMES]


if __name__ == "__main__":
    run()
