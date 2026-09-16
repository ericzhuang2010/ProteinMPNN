"""
aibinder_tools.py — shared helpers for the ImmunoVerse agent tools.

NOTE: In the real ImmunoVerse repo this module is much larger (it holds the
analyze_pdb_interactions tool, etc.). This is a minimal, standalone stand-in
that provides only the shared conventions other tool modules reuse — most
importantly _find_pdb — so proteinmpnn_tools.py can do:

    from aibinder_tools import _find_pdb

Conventions modeled here (the ones the assignment points at):
  * Paths come from env vars with sensible defaults (BASE_DIR).
  * _find_pdb(pdb_path) resolves a PDB so callers can pass just a filename.
  * Heavy libs are imported lazily so importing this module never crashes.
"""

import os
from pathlib import Path
from typing import Any, List, Optional

# ---- paths from env vars with defaults --------------------------------------
# Root under which bare PDB filenames are resolved. Override with AIBINDER_BASE_DIR.
BASE_DIR = os.environ.get("AIBINDER_BASE_DIR", os.getcwd())

# Subdirectories of BASE_DIR commonly searched for structure files.
_PDB_SUBDIRS = ("", "inputs", "pdbs", "structures", "data")

# ---- lazy heavy imports (never crash on import) -----------------------------
_NUMPY_OK = False
try:
    import numpy as _np  # noqa: F401

    _NUMPY_OK = True
except ImportError:
    _np = None


def _find_pdb(pdb_path: str) -> Optional[str]:
    """Resolve a PDB path so callers can pass just a filename.

    Resolution order:
      1. The path exactly as given (absolute or relative to the CWD).
      2. The same, with a ``.pdb`` extension appended if none was provided.
      3. Each candidate name under BASE_DIR and a few common subdirectories
         (inputs/, pdbs/, structures/, data/).

    Returns an absolute path to an existing file, or None if nothing matches.
    """
    if not pdb_path:
        return None

    given = Path(pdb_path).expanduser()
    # Candidate filenames: as given, plus a ".pdb"-suffixed variant if bare.
    candidate_names: List[str] = [str(given)]
    if not given.suffix:
        candidate_names.append(str(given) + ".pdb")

    search_roots = [Path.cwd(), Path(BASE_DIR)]
    for name in candidate_names:
        cand = Path(name).expanduser()
        if cand.is_file():
            return str(cand.resolve())
        # Only try BASE_DIR subdirs for non-absolute names.
        if not cand.is_absolute():
            for root in search_roots:
                for sub in _PDB_SUBDIRS:
                    probe = (root / sub / name) if sub else (root / name)
                    if probe.is_file():
                        return str(probe.resolve())
    return None

