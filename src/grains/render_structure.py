"""Render a gluten peptide sitting in the HLA-DQ2.5 (or DQ8) binding groove.

CPU-only, no external viewer required (no ChimeraX / PyMOL install, which this machine
can't spare the disk for): parses the PDB file with Biopython and draws a CA backbone
trace with matplotlib's 3D axes. Not a publication-grade cartoon rendering, but enough
to show the peptide nestled between the two MHC chains for Checkpoint 1.

Usage:
    PYTHONPATH=src python -m grains.render_structure 1S9V A B C \
        "HLA-DQ2.5 + DQ2.5-glia-a1a (1S9V)" results/figures/1S9V_dq25_groove.png
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from Bio.PDB import PDBParser  # noqa: E402

RAW_PDB_DIR = Path("data/raw/pdb")

CHAIN_STYLE = {
    "alpha": dict(color="#8fa3b3", label="MHC alpha chain", lw=1.6),
    "beta": dict(color="#4c6a82", label="MHC beta chain", lw=1.6),
    "peptide": dict(color="#d1495b", label="peptide", lw=4.0),
}


def ca_trace(structure, chain_id: str):
    model = next(structure.get_models())
    chain = model[chain_id]
    coords = [
        residue["CA"].get_coord()
        for residue in chain
        if residue.has_id("CA")
    ]
    xs, ys, zs = zip(*coords)
    return xs, ys, zs


def render(pdb_id: str, alpha_chain: str, beta_chain: str, peptide_chain: str,
           title: str, out_path: Path) -> None:
    pdb_file = RAW_PDB_DIR / f"{pdb_id}.pdb"
    structure = PDBParser(QUIET=True).get_structure(pdb_id, pdb_file)

    fig = plt.figure(figsize=(7, 6))
    ax = fig.add_subplot(projection="3d")

    for chain_id, role in [(alpha_chain, "alpha"), (beta_chain, "beta"), (peptide_chain, "peptide")]:
        xs, ys, zs = ca_trace(structure, chain_id)
        style = CHAIN_STYLE[role]
        ax.plot(xs, ys, zs, color=style["color"], linewidth=style["lw"], label=style["label"])
        if role == "peptide":
            ax.scatter(xs, ys, zs, color=style["color"], s=18)

    ax.set_title(f"{title}\nCA backbone trace, chains {alpha_chain}/{beta_chain}/{peptide_chain}")
    ax.set_xlabel("x (Å)")
    ax.set_ylabel("y (Å)")
    ax.set_zlabel("z (Å)")
    ax.legend(loc="upper left")
    ax.view_init(elev=20, azim=-60)
    fig.tight_layout()

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=200)
    print(f"wrote {out_path}")


def main(argv: list[str]) -> None:
    pdb_id, alpha_chain, beta_chain, peptide_chain, title, out = argv[1:7]
    render(pdb_id, alpha_chain, beta_chain, peptide_chain, title, Path(out))


if __name__ == "__main__":
    main(sys.argv)
