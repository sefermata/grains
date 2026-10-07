"""Build PANDORA's MHC class II template database (Checkpoint 3 prerequisite).

PANDORA models a peptide-MHC complex by homology: it finds experimentally solved complexes of
the same MHC class, picks the best template, and threads the target peptide onto it. So before
anything can be modeled, that library of solved structures has to exist locally. This module
builds it.

Why this isn't just `Database.construct_database(download=True)`
----------------------------------------------------------------
PANDORA 0.9's own downloader cannot run here at all:

1. It requests IMGT over `http://`, and `www.imgt.org` is https-only now (plain HTTP returns
   nothing — not a redirect, so there is nothing to follow).
2. It shells out to `wget`, which macOS does not ship.
3. It does `gunzip archive.tgz` then `tar -xvf archive.tar`, which leaves a fully uncompressed
   tar sitting beside the extracted tree. The archive is 1.75 GB compressed; that peak does not
   fit this machine's free space.

So the acquisition is done here instead, and `construct_database` is called with
`download=False` to do only the parsing.

Two deliberate reductions, both large:

- **Class II only.** This project models HLA-DQ2.5, which is class II. Parsing class I as well
  would roughly double the time and disk for templates that can never be selected.
- **Only the class II members are extracted from the archive.** IMGT's flat-file set covers
  every immune structure it holds — antibodies, TCRs, MHC class I — as one `IMGT-<id>.pdb.gz`
  per entry. Asking IMGT which PDB IDs are class II *first* (205 of them) means extracting 205
  members instead of thousands. The 1.75 GB download is unavoidable (a `.tar.gz` can't be
  randomly accessed), but it is streamed, never expanded in full, and deleted afterward.

Usage
-----
    PYTHONPATH=src python -m grains.pandora_db           # build (resumable)
    PYTHONPATH=src python -m grains.pandora_db --keep    # keep the archive for a re-run

Output
------
    <PANDORA_data>/mhcII_database.pkl    the pickled template database
    data/processed/pandora_db_build.json  time, disk, and template-count accounting
"""

from __future__ import annotations

import json
import shutil
import subprocess
import sys
import tarfile
import time
from pathlib import Path

from . import config, pandora_compat  # noqa: F401  (compat must precede PANDORA)

import PANDORA  # noqa: E402

IMGT_ARCHIVE_URL = "https://www.imgt.org/download/3Dstructure-DB/IMGT3DFlatFiles.tgz"
MEMBER_DIR = "IMGT3DFlatFiles"

DATA_DIR = Path(PANDORA.PANDORA_data)
DATABASE_PATH = DATA_DIR / "mhcII_database.pkl"
ARCHIVE_PATH = DATA_DIR / "PDBs" / "IMGT_retrieved" / "IMGT3DFlatFiles.tgz"
REPORT_PATH = config.PROCESSED / "pandora_db_build.json"

# The tree PANDORA expects to already exist. It ships none of it, and its functions `os.chdir`
# into these paths without creating them, so a missing directory is a FileNotFoundError.
REQUIRED_DIRS = (
    "PDBs/IMGT_retrieved",
    "PDBs/pMHCI",
    "PDBs/pMHCII",
    "PDBs/Bad/pMHCI",
    "PDBs/Bad/pMHCII",
    "csv_pkl_files",
    "outputs",          # Pandora.__log writes Pandora_log.txt here after a successful model
)


def disk_free_gb() -> float:
    return shutil.disk_usage("/").free / 1e9


def dir_size_gb(path: Path) -> float:
    if not path.exists():
        return 0.0
    return sum(f.stat().st_size for f in path.rglob("*") if f.is_file()) / 1e9


def ensure_tree(data_dir: Path = DATA_DIR) -> None:
    """Create the directory layout PANDORA assumes but does not ship."""
    for relative in REQUIRED_DIRS:
        (data_dir / relative).mkdir(parents=True, exist_ok=True)


def class_ii_ids() -> list[str]:
    """PDB IDs of every MHC class II structure IMGT holds, over https."""
    return pandora_compat.download_ids_imgt("MH2")


def download_archive(dest: Path = ARCHIVE_PATH, url: str = IMGT_ARCHIVE_URL) -> Path:
    """Fetch the IMGT flat-file archive with curl, resuming a partial download.

    curl rather than urllib so a 1.75 GB transfer is resumable across interruptions, and
    rather than wget because macOS ships curl and not wget.
    """
    dest.parent.mkdir(parents=True, exist_ok=True)
    print(f"downloading {url}", flush=True)
    print(f"  -> {dest}", flush=True)
    subprocess.run(
        ["curl", "-fL", "--retry", "3", "--retry-delay", "5", "-C", "-",
         "-o", str(dest), url],
        check=True,
    )
    print(f"  archive {dest.stat().st_size / 1e9:.2f} GB", flush=True)
    return dest


def extract_members(archive: Path, pdb_ids: list[str], data_dir: Path = DATA_DIR) -> tuple[int, list[str]]:
    """Extract only `IMGT-<id>.pdb.gz` for the given IDs. Returns (written, missing).

    Streams the archive once. Members are matched against a set, so this is one pass however
    many IDs are requested, and the uncompressed tar is never written to disk.
    """
    out_dir = data_dir / "PDBs" / "IMGT_retrieved" / MEMBER_DIR
    out_dir.mkdir(parents=True, exist_ok=True)

    wanted = {f"IMGT-{pdb_id}.pdb.gz": pdb_id for pdb_id in pdb_ids}
    found: set[str] = set()
    written = 0

    with tarfile.open(archive, "r|gz") as tar:          # "r|gz" = stream, no seeking
        for member in tar:
            if not member.isfile():
                continue
            name = Path(member.name).name
            if name not in wanted:
                continue
            source = tar.extractfile(member)
            if source is None:
                continue
            with open(out_dir / name, "wb") as target:
                shutil.copyfileobj(source, target)
            found.add(name)
            written += 1
            if written % 25 == 0:
                print(f"  extracted {written}/{len(wanted)}", flush=True)

    missing = sorted(wanted[name] for name in wanted.keys() - found)
    print(f"  extracted {written}/{len(wanted)} class II structures", flush=True)
    if missing:
        # Worth stating rather than swallowing: IMGT's query and its archive are published
        # separately and can disagree, and a template silently absent from the library is
        # exactly the kind of gap that shows up later as an unexplained modeling failure.
        print(f"  WARNING: {len(missing)} ID(s) in IMGT's class II list are absent from the "
              f"archive: {', '.join(missing)}", flush=True)
    return written, missing


def build(keep_archive: bool = False) -> dict:
    """Acquire the class II templates and construct the database. Returns the build report."""
    pandora_compat.install_imgt_https()
    from PANDORA.Database import Database

    start = time.time()
    free_before = disk_free_gb()
    print(f"data_dir   {DATA_DIR}")
    print(f"disk free  {free_before:.1f} GB\n")

    ensure_tree()

    pdb_ids = class_ii_ids()
    print(f"IMGT lists {len(pdb_ids)} MHC class II structures\n")

    if not ARCHIVE_PATH.exists():
        download_archive()
    else:
        print(f"archive already present ({ARCHIVE_PATH.stat().st_size / 1e9:.2f} GB), "
              "skipping download")

    print("\nextracting class II members...", flush=True)
    written, missing = extract_members(ARCHIVE_PATH, pdb_ids)

    if not keep_archive:
        ARCHIVE_PATH.unlink(missing_ok=True)
        print(f"removed the archive; disk free {disk_free_gb():.1f} GB")

    print("\nparsing class II structures into a template database...", flush=True)
    database = Database.Database()
    database.construct_database(
        save=str(DATABASE_PATH),
        data_dir=str(DATA_DIR),
        MHCI=False,
        MHCII=True,
        download=False,          # acquisition already done above, over https
        update_ref_sequences=False,  # class II targets must supply their own chain sequences
    )

    elapsed = time.time() - start
    report = {
        "built_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "imgt_class_ii_listed": len(pdb_ids),
        "members_extracted": written,
        "members_missing_from_archive": missing,
        "templates_parsed": len(database.MHCII_data),
        "wall_clock_min": round(elapsed / 60, 1),
        "data_dir_gb": round(dir_size_gb(DATA_DIR), 2),
        "disk_free_before_gb": round(free_before, 1),
        "disk_free_after_gb": round(disk_free_gb(), 1),
        "database_path": str(DATABASE_PATH),
    }

    config.ensure_dirs()
    REPORT_PATH.write_text(json.dumps(report, indent=2) + "\n")

    print("\n=== database build done ===")
    print(f"templates   {report['templates_parsed']} of {written} extracted")
    print(f"wall clock  {report['wall_clock_min']} min")
    print(f"data_dir    {report['data_dir_gb']} GB")
    print(f"disk free   {free_before:.1f} -> {report['disk_free_after_gb']} GB")
    print(f"saved       {DATABASE_PATH}")
    print(f"report      {REPORT_PATH}")
    return report


#: IMGT's flat files, where the TCR annotation survives. PANDORA's cleaned templates keep only
#: chains M/N/P, so by the time a structure is a Template object you can no longer tell whether
#: it was crystallized with a T-cell receptor bound.
IMGT_FLAT_DIR = DATA_DIR / "PDBs" / "IMGT_retrieved" / MEMBER_DIR


def is_tcr_bound(pdb_id: str, flat_dir: Path = IMGT_FLAT_DIR) -> bool | None:
    """Whether this structure was crystallized with a T-cell receptor bound (ternary).

    Read from IMGT's own chain annotation: `REMARK 410 IMGT chain description  TR-ALPHA`
    and friends. Returns None when the flat file isn't on disk.

    This matters because a homology model inherits its template's peptide conformation, and
    the peptide sits differently in a TCR-bound complex than in a bare peptide-MHC one.
    """
    import gzip

    path = flat_dir / f"IMGT-{pdb_id}.pdb.gz"
    if not path.exists():
        return None
    with gzip.open(path, "rt", errors="replace") as handle:
        for line in handle:
            if not line.startswith("REMARK 410"):
                if line.startswith(("ATOM", "HETATM")):
                    break          # past the remarks; stop reading the coordinates
                continue
            if "IMGT chain description" in line:
                description = line.split("IMGT chain description", 1)[1].strip()
                if description.startswith("TR-"):
                    return True
    return False


def tcr_bound(database, flat_dir: Path = IMGT_FLAT_DIR) -> set[str]:
    """IDs of every template in `database` crystallized with a TCR bound."""
    return {pdb_id for pdb_id in database.MHCII_data
            if is_tcr_bound(pdb_id, flat_dir) is True}


def database_ready() -> tuple[bool, str]:
    """Whether a built class II template database is on disk."""
    if not DATABASE_PATH.exists():
        return False, (f"no template database at {DATABASE_PATH}. Build it with:\n"
                       "    PYTHONPATH=src python -m grains.pandora_db")
    size_mb = DATABASE_PATH.stat().st_size / 1e6
    return True, f"{DATABASE_PATH.name} present ({size_mb:.0f} MB)"


def main() -> None:
    build(keep_archive="--keep" in sys.argv)


if __name__ == "__main__":
    main()
