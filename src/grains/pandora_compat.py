"""Make PANDORA 0.9 importable against Biopython 1.88.

PANDORA (`csb-pandora` 0.9, internal version 1.0) was written against Biopython < 1.80 and
calls three APIs that have since been removed. 0.9 is the only release on PyPI, so there is no
upstream fix to upgrade to:

| PANDORA calls | Removed in | Used by |
|---|---|---|
| `Bio.PDB.Polypeptide.three_to_one` | 1.80 | `PMHC/Model.py` `renumber()` |
| `Bio.pairwise2` | 1.84 | `renumber()`, and the NetMHCpan anchor path |
| `Bio.Align.Applications.MuscleCommandline` | 1.82 | `Pandora/Align.py`, all MSA calls |

Without these three names the import fails at module load, so nothing runs at all — not even
the database build.

**Separately, PANDORA's data acquisition is also dead.** It requests IMGT over `http://`, and
`www.imgt.org` is now https-only — plain HTTP returns nothing at all (curl exit 000, not a
redirect, so there is nothing for urllib to follow). That breaks both places it reaches out:

| PANDORA calls | Problem |
|---|---|
| `download_unzip_imgt_structures` | `http://` URL, and shells out to `wget` (absent on macOS) |
| `download_ids_imgt` | `http://` URL — and this one runs even with `download=False` |

`download_ids_imgt` is replaced here with an https equivalent, parsing the response exactly as
the original did so the returned ID list is identical. The structure archive is handled by
`grains.pandora_db` instead, which also avoids PANDORA's `gunzip`-then-`tar` sequence — that
materializes a full uncompressed tar *beside* the extracted tree, and at 1.75 GB compressed the
peak does not fit this machine's free space.

**Import this module before importing PANDORA.** `install()` is idempotent and runs on import.

Design decision: the two shims that are *off this project's code path* raise
`NotImplementedError` when actually called, rather than being reimplemented.

`three_to_one` and `pairwise2` are reachable only from `Model.renumber()` (a benchmarking
helper) and from NetMHCpan-based anchor prediction. This project passes `benchmark=False` and
supplies anchors explicitly, so neither runs. Reimplementing them would mean shipping an
approximation of an alignment into a structural modeling pipeline on a path that is never
exercised or tested — precisely the "silently plausible wrong number" this project keeps
paying for. A loud crash on an unexpected path is the safer failure.

`MuscleCommandline` is different: it *is* on the class II path (`Align.py` aligns the M and N
MHC chains against the templates with it), so it has to genuinely work. PANDORA uses it as

    msl = MuscleCommandline(input=..., out=...)
    os.system(str(msl) + ' -quiet')

i.e. it only ever needs `str(msl)` to be a runnable command line. The shim reproduces
Biopython's MUSCLE v3 syntax (`muscle -in IN -out OUT`), which is why **MUSCLE 3.8.1551** is
the required version — v5 renamed those flags to `-align`/`-output` and would silently produce
no output file here.

    conda install -c bioconda muscle=3.8.1551
"""

from __future__ import annotations

import shlex
import shutil
import sys
import types

MUSCLE_BINARY = "muscle"

#: The MUSCLE major version whose command-line syntax `MuscleShim` emits.
MUSCLE_SYNTAX_VERSION = 3


def _three_to_one(resname: str) -> str:
    """Removed in Biopython 1.80. Off this project's path -- see the module docstring."""
    raise NotImplementedError(
        f"three_to_one({resname!r}) was called, but grains.pandora_compat only stubs it.\n"
        "It is reachable from PANDORA's Model.renumber() and NetMHCpan anchor prediction, "
        "neither of which this project uses (benchmark=False, anchors supplied explicitly).\n"
        "Reaching it means the code path changed. Implement it with "
        "Bio.Data.IUPACData.protein_letters_3to1 -- and add a test that covers the new path."
    )


def _pairwise2_globalxx(*args, **kwargs):
    """Removed in Biopython 1.84. Off this project's path -- see the module docstring."""
    raise NotImplementedError(
        "pairwise2.align.globalxx was called, but grains.pandora_compat only stubs it.\n"
        "It is reachable from PANDORA's Model.renumber() and NetMHCpan anchor prediction, "
        "neither of which this project uses (benchmark=False, anchors supplied explicitly).\n"
        "Reaching it means the code path changed. Bio.Align.PairwiseAligner reproduces "
        "globalxx exactly with match=1, mismatch=0, open_gap_score=0, extend_gap_score=0, "
        "mode='global' -- and add a test that covers the new path."
    )


class MuscleShim:
    """Stand-in for `Bio.Align.Applications.MuscleCommandline`, MUSCLE v3 syntax.

    PANDORA only ever renders this to a string and hands it to `os.system`, so `__str__` is
    the whole contract. Paths are shell-quoted; PANDORA appends ` -quiet` to the result, which
    stays outside the quoting.
    """

    def __init__(self, cmd: str = MUSCLE_BINARY, *, input: str = "", out: str = "", **kwargs):
        self.cmd = cmd
        self.input = input
        self.out = out
        self.extra = kwargs

    def __str__(self) -> str:
        executable = shutil.which(self.cmd)
        if executable is None:
            raise FileNotFoundError(
                f"MUSCLE ({self.cmd!r}) is not on PATH, and PANDORA needs it to align the "
                "MHC class II M and N chains against its templates.\n"
                "    conda install -c bioconda muscle=3.8.1551\n"
                f"Version {MUSCLE_SYNTAX_VERSION}.x specifically: v5 renamed -in/-out to "
                "-align/-output and would write no output file while still exiting 0."
            )
        return f"{shlex.quote(executable)} -in {shlex.quote(self.input)} -out {shlex.quote(self.out)}"

    def __repr__(self) -> str:
        return f"MuscleShim(input={self.input!r}, out={self.out!r})"


def muscle_version() -> str | None:
    """The installed MUSCLE's version banner, or None if it isn't on PATH."""
    import subprocess

    executable = shutil.which(MUSCLE_BINARY)
    if executable is None:
        return None
    try:
        result = subprocess.run([executable, "-version"], capture_output=True, text=True,
                                timeout=30)
    except (OSError, subprocess.SubprocessError):
        return None
    return (result.stdout or result.stderr).strip().splitlines()[0] if (
        result.stdout or result.stderr) else None


def muscle_syntax_ok() -> tuple[bool, str]:
    """Whether the installed MUSCLE speaks the `-in`/`-out` syntax the shim emits."""
    banner = muscle_version()
    if banner is None:
        return False, ("MUSCLE is not installed (conda install -c bioconda muscle=3.8.1551). "
                       "PANDORA cannot align MHC class II chains without it.")
    if "v3." not in banner and "3.8" not in banner:
        return False, (
            f"MUSCLE reports {banner!r}, but the shim emits v3 `-in`/`-out` syntax. "
            "v5 renamed these to `-align`/`-output`, so alignments would silently not be "
            "written. Install muscle=3.8.1551."
        )
    return True, banner


IMGT_QUERY_URL = "https://www.imgt.org/3Dstructure-DB/cgi/3Dquery.cgi"


def download_ids_imgt(ReceptorType, data_dir=None, out_tsv=False):
    """https replacement for `PANDORA.Database.Database_functions.download_ids_imgt`.

    Same arguments, same return value, same optional TSV side effect. The response is parsed
    exactly as the original parsed it, so the ID list is identical to what PANDORA would have
    produced back when its http URL still worked.
    """
    import urllib.parse
    import urllib.request

    import PANDORA

    if data_dir is None:
        data_dir = PANDORA.PANDORA_data

    payload = urllib.parse.urlencode(
        {"ReceptorType": ReceptorType, "type-entry": "PDB"}
    ).encode("ascii")
    request = urllib.request.Request(IMGT_QUERY_URL, payload)
    with urllib.request.urlopen(request, timeout=180) as response:
        text = response.read().decode("utf-8").splitlines()

    # Parsing kept byte-for-byte identical to PANDORA 0.9's own.
    ids = [x for x in text if "href" in x and "pdbcode" in x]
    ids = [x.split('"') for x in ids]
    ids = [x[3][-4:] for x in ids]

    if out_tsv:
        out_path = f"{data_dir}/csv_pkl_files/{out_tsv}"
        with open(out_path, "w") as handle:
            handle.write(ReceptorType + " IMGT IDs\n")
            for identifier in ids:
                handle.write(identifier + "\n")

    return ids


def install() -> None:
    """Inject the three removed names into Biopython. Idempotent."""
    from Bio.PDB import Polypeptide

    if not hasattr(Polypeptide, "three_to_one"):
        Polypeptide.three_to_one = _three_to_one

    if "Bio.pairwise2" not in sys.modules:
        pairwise2 = types.ModuleType("Bio.pairwise2")
        align = types.SimpleNamespace(
            globalxx=_pairwise2_globalxx,
            globalds=_pairwise2_globalxx,
            globalms=_pairwise2_globalxx,
        )
        pairwise2.align = align
        sys.modules["Bio.pairwise2"] = pairwise2
        import Bio

        Bio.pairwise2 = pairwise2

    if "Bio.Align.Applications" not in sys.modules:
        applications = types.ModuleType("Bio.Align.Applications")
        applications.MuscleCommandline = MuscleShim
        sys.modules["Bio.Align.Applications"] = applications
        import Bio.Align

        Bio.Align.Applications = applications


def install_imgt_https() -> None:
    """Point PANDORA's IMGT ID query at https. Separate from `install()` because it has to
    run *after* PANDORA is importable, which `install()` is what makes possible."""
    from PANDORA.Database import Database, Database_functions

    Database_functions.download_ids_imgt = download_ids_imgt
    # Database.py did `from PANDORA.Database import Database_functions`, so it resolves the
    # attribute at call time off that module object -- patching the module is enough. Assert
    # it rather than trusting it: a future PANDORA doing `from ... import download_ids_imgt`
    # would bind the old function directly and silently keep using http.
    assert Database.Database_functions.download_ids_imgt is download_ids_imgt, (
        "patching Database_functions.download_ids_imgt did not take effect; PANDORA's "
        "import style must have changed"
    )


install()
