"""Minimal shared PDB-path helper used by the ProteinMPNN tool delivery.

In ImmunoVerse this helper comes from the full ``aibinder_tools`` module.  A
small copy is included here so the milestone delivery and its tests can run on
their own.
"""

import os
from pathlib import Path
from typing import List, Optional


BASE_DIR = os.environ.get("AIBINDER_BASE_DIR", os.getcwd())
_PDB_SUBDIRS = ("", "inputs", "pdbs", "structures", "data")


def _find_pdb(pdb_path: str) -> Optional[str]:
    """Return the absolute path to a PDB, or ``None`` when it cannot be found."""
    if not pdb_path:
        return None

    given = Path(pdb_path).expanduser()
    candidate_names: List[str] = [str(given)]
    if not given.suffix:
        candidate_names.append(str(given) + ".pdb")

    search_roots = [Path.cwd(), Path(BASE_DIR)]
    for name in candidate_names:
        candidate = Path(name).expanduser()
        if candidate.is_file():
            return str(candidate.resolve())
        if not candidate.is_absolute():
            for root in search_roots:
                for subdir in _PDB_SUBDIRS:
                    probe = root / subdir / name if subdir else root / name
                    if probe.is_file():
                        return str(probe.resolve())
    return None
