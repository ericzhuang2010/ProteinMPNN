"""
proteinmpnn_tools.py — Agent-callable wrappers for ProteinMPNN
(protein sequence design & scoring from a 3D backbone).
"""
import os, re, json, tempfile, subprocess, shutil
from pathlib import Path
from typing import Any, Dict, List

# ---- locate ProteinMPNN (no hardcoding) -------------------------------
PROTEINMPNN_DIR = os.environ.get("PROTEINMPNN_DIR", "/opt/ProteinMPNN")
RUN_SCRIPT = os.path.join(PROTEINMPNN_DIR, "protein_mpnn_run.py")

# ---- lazy dependency check (copy the aibinder_tools pattern) ----------
_TORCH_OK = False
try:
    import torch  # noqa
    _TORCH_OK = True
except ImportError:
    pass

# numpy is needed to read score-only .npz files (it ships with torch anyway).
_NUMPY_OK = False
try:
    import numpy as np  # noqa
    _NUMPY_OK = True
except ImportError:
    np = None


def _check_deps() -> Dict[str, Any]:
    """Return an error dict if something is missing, else None."""
    missing = []
    if not _TORCH_OK:
        missing.append("torch")
    if not _NUMPY_OK:
        missing.append("numpy")
    if not os.path.isfile(RUN_SCRIPT):
        missing.append(f"ProteinMPNN (set PROTEINMPNN_DIR; looked in {PROTEINMPNN_DIR})")
    if missing:
        return {"success": False,
                "error": "Missing dependencies: " + ", ".join(missing) +
                ". Install torch and clone ProteinMPNN."}
    return None


# reuse the existing PDB resolver so users can pass just a filename
from aibinder_tools import _find_pdb  # noqa


# =======================================================================
# Parser — write this FIRST and unit-test it (no torch needed)
# =======================================================================
# Different types of header fields.
_FLOAT_KEYS = {"score", "global_score", "seq_recovery", "T"}
_INT_KEYS = {"sample", "seed"}
# Split header fields on commas that are NOT inside [bracketed, lists]: chain
# lists like fixed_chains=['A', 'B', 'C'] contain commas that are part of the
# value, not field separators. The negative lookahead "(?![^\[]*\])" means
# "this comma is not followed by ...] before the next '[' ", i.e. not inside [].
_FIELD_SPLIT_RE = re.compile(r",(?![^\[]*\])")


def _parse_header_value(key: str, value: str) -> Any:
    """Convert a raw header value to float/int/list based on its key."""
    value = value.strip()
    if value.startswith("["):  # e.g. "['A', 'B', 'C']" -> ["A", "B", "C"]
        return [p.strip().strip("'\"") for p in value.strip("[]").split(",") if p.strip()]
    if key in _FLOAT_KEYS:
        try:
            return float(value)
        except ValueError:
            # If the value is not a number, return it as is ?
            return value
    if key in _INT_KEYS:
        try:
            return int(value)
        except ValueError:
            # If the value is not a number, return it as is ?
            return value
    return value


def _parse_mpnn_fasta(text: str) -> List[Dict[str, Any]]:
    """
    Parse ProteinMPNN's output FASTA into a list of records:
    [{"header": {...}, "sequence": "..."}, ...]
    The first record is the native sequence. See assignment Appendix C.
    HINT: split into >header / sequence pairs, then split the header on
    commas, then each piece on the first '=' to get key=value pairs.
    """
    records = []
    blocks = [b for b in text.strip().split(">") if b.strip()]
    for b in blocks:
        lines = b.splitlines()
        header_line = lines[0]
        seq = "".join(lines[1:]).strip()
        meta: Dict[str, Any] = {}
        for i, piece in enumerate(_FIELD_SPLIT_RE.split(header_line)):
            piece = piece.strip()
            if "=" in piece:
                # key=value field; values are coerced to float/int/list.
                k, v = piece.split("=", 1)
                meta[k.strip()] = _parse_header_value(k.strip(), v)
            elif i == 0 and piece:
                # The first record's name has no '=' before the first comma.
                meta["name"] = piece
        records.append({"header": meta, "sequence": seq})
    return records


# =======================================================================
# Milestone 2 — design
# =======================================================================
def proteinmpnn_design_sequences(pdb_path: str, chains_to_design: str = "",
        num_sequences: int = 8, sampling_temp: float = 0.1,
        seed: int = 37, model_name: str = "v_48_020") -> Dict[str, Any]:
    """Design amino-acid sequences for a protein backbone using ProteinMPNN.

    Designs the requested chain(s) while keeping the rest as fixed context, and
    returns the native sequence + score and a list of designs with score,
    global_score, seq_recovery, sample and temperature. Score is a negative
    log-likelihood; lower = more compatible. Returns the Appendix-A dict, or
    {"success": False, "error": "..."} (never raises).
    """
    dep = _check_deps()
    if dep: return dep
    resolved = _find_pdb(pdb_path)
    if not resolved:
        return {"success": False, "error": f"PDB file not found: {pdb_path}"}
    out_dir = tempfile.mkdtemp(prefix="mpnn_")
    try:
        cmd = ["python", RUN_SCRIPT,
               "--pdb_path", resolved,
               "--out_folder", out_dir,
               "--num_seq_per_target", str(num_sequences),
               "--sampling_temp", str(sampling_temp),
               "--seed", str(seed),
               "--model_name", model_name,
               # batch_size=1: CPU-friendly; batching only speeds up GPU runs.
               "--batch_size", "1"]
        if chains_to_design.strip():
            cmd += ["--pdb_path_chains", chains_to_design.strip()]
        
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
        if proc.returncode != 0:
            return {"success": False,
                    "error": f"ProteinMPNN failed: {proc.stderr[-800:]}"}
        stem = Path(resolved).stem
        fa = Path(out_dir) / "seqs" / f"{stem}.fa"
        if not fa.is_file():
            return {"success": False,
                    "error": f"ProteinMPNN produced no FASTA at {fa}. {proc.stderr[-800:]}"}
        records = _parse_mpnn_fasta(fa.read_text())
        if not records:
            return {"success": False, "error": f"Could not parse {fa}."}

        # Keep the .fa for debugging (like AIbinder keeps the .pml): copy it to a
        # more stable path, then clean up the temp folder in the finally block.
        keep_fa = tempfile.NamedTemporaryFile(prefix=f"{stem}_mpnn_", suffix=".fa", delete=False)
        keep_fa.close()
        shutil.copyfile(fa, keep_fa.name)

        # Build the result dict EXACTLY as in assignment Appendix A.
        native = records[0]
        native_header = native["header"]
        sequences = [{
            "sample": d["header"].get("sample"),
            "sequence": d["sequence"],
            "score": d["header"].get("score"),
            "global_score": d["header"].get("global_score"),
            "seq_recovery": d["header"].get("seq_recovery"),
            "temperature": d["header"].get("T"),
        } for d in records[1:]]
        return {
            "success": True,
            "pdb_file": resolved,
            "designed_chains": native_header.get("designed_chains", []),
            "fixed_chains": native_header.get("fixed_chains", []),
            "model_name": native_header.get("model_name", model_name),
            "sampling_temp": float(sampling_temp),
            "num_sequences": int(num_sequences),
            "native_sequence": native["sequence"],
            "native_score": native_header.get("score"),
            "sequences": sequences,
            "out_fasta_path": keep_fa.name,
        }
    except subprocess.TimeoutExpired:
        return {"success": False, "error": "ProteinMPNN timed out (300s)."}
    except Exception as exc:  # never raise out of a public tool function
        return {"success": False, "error": f"ProteinMPNN design error: {exc}"}
    finally:
        shutil.rmtree(out_dir, ignore_errors=True)


# =======================================================================
# Milestone 3 — score-only (implemented with --score_only 1)
# =======================================================================
# Where the score is read from: --score_only both prints a rounded summary to
# stdout AND writes .npz files under <out>/score_only/. We read the .npz files
# because they expose the raw, full-precision `score`/`global_score` arrays
# (stdout is rounded to 4 dp), and reading named arrays from a file is more
# robust than regex-parsing log lines. See docs/score_output_choice.md.


def _read_fasta_names(fasta_file: str) -> List[str]:
    """Return record names (text after '>') from a FASTA file, in order."""
    names: List[str] = []
    try:
        for line in Path(fasta_file).read_text().splitlines():
            if line.startswith(">"):
                names.append(line[1:].strip() or f"seq{len(names) + 1}")
    except OSError:
        pass
    return names


def _score_from_npz(npz_path: Path, name: str) -> Dict[str, Any]:
    """Load one score-only .npz; return {name, score, global_score}.

    Averages the per-pass arrays and casts to float so the result stays
    JSON-serializable (no numpy scalars leak out).
    """
    data = np.load(str(npz_path))
    return {
        "name": name,
        "score": float(data["score"].mean()),
        "global_score": float(data["global_score"].mean()),
    }


def _collect_scores(score_dir: Path, stem: str,
                    fasta_seq_names: List[str]) -> List[Dict[str, Any]]:
    """Read score-only .npz files into [{name, score, global_score}].

    score_only writes <stem>_pdb.npz (the backbone's native sequence) and
    <stem>_fasta_<fc>.npz (fc = 1, 2, ... for each input FASTA sequence). The
    fc index maps back to fasta_seq_names, the names of the sequences we scored.
    """
    scored: List[Dict[str, Any]] = []
    pdb_npz = score_dir / f"{stem}_pdb.npz"
    if pdb_npz.is_file():
        scored.append(_score_from_npz(pdb_npz, f"{stem} (native from PDB)"))
    fc = 1
    while True:
        fasta_npz = score_dir / f"{stem}_fasta_{fc}.npz"
        if not fasta_npz.is_file():
            break
        if fc - 1 < len(fasta_seq_names):
            name = fasta_seq_names[fc - 1]
        else:
            name = f"fasta_{fc}"
        scored.append(_score_from_npz(fasta_npz, name))
        fc += 1
    return scored


def proteinmpnn_score_sequence(pdb_path: str, sequence: str = "",
        fasta_path: str = "", chains_to_score: str = "", seed: int = 37) -> Dict[str, Any]:
    """Score sequence(s) against a PDB backbone using ProteinMPNN (score-only).

    Provide a single ``sequence`` string, a ``fasta_path`` of sequences, or
    neither (scores only the backbone's own native sequence). Returns the
    Appendix-A dict, or {"success": False, "error": "..."} (never raises).
    Score is a negative log-likelihood; lower means more compatible. The score
    is read from ProteinMPNN's score_only/*.npz files (full precision).
    """
    dep = _check_deps()
    if dep: return dep
    resolved = _find_pdb(pdb_path)
    if not resolved:
        return {"success": False, "error": f"PDB file not found: {pdb_path}"}
    out_dir = tempfile.mkdtemp(prefix="mpnn_")
    try:
        # Normalize the caller's input to a single --path_to_fasta value. If both
        # fasta_path and sequence are given, use fasta_path. If neither is
        # given, path_to_fasta stays "" and only the PDB's native sequence is scored.
        # fasta_seq_names tracks the input sequence names in order, so the
        # score_only/<stem>_fasta_<fc>.npz files can be mapped back to a name.
        path_to_fasta = ""
        fasta_seq_names: List[str] = []
        if fasta_path:
            resolved_fasta = fasta_path if os.path.isfile(fasta_path) else _find_pdb(fasta_path)
            if not resolved_fasta or not os.path.isfile(resolved_fasta):
                return {"success": False, "error": f"FASTA file not found: {fasta_path}"}
            path_to_fasta = os.path.abspath(resolved_fasta)
            fasta_seq_names = _read_fasta_names(path_to_fasta)
        elif sequence.strip():
            tmp_fa = Path(out_dir) / "input_to_score.fa"
            tmp_fa.write_text(f">seq1\n{sequence.strip()}\n")
            path_to_fasta = str(tmp_fa)
            fasta_seq_names = ["seq1"]

        cmd = ["python", RUN_SCRIPT,
               "--pdb_path", resolved,
               "--out_folder", out_dir,
               "--score_only", "1",
               "--seed", str(seed),
               # num_seq_per_target=1: in score-only this is NUM_BATCHES, the
               # number of random decoding-order passes averaged into the score.
               # 1 = a single scoring pass per sequence (set >1 to average for a
               # more stable mean). Pinned explicitly so we don't depend on the
               # run script's default value.
               "--num_seq_per_target", "1",
               # batch_size=1: CPU-friendly (batching only helps GPU), and with
               # num_seq_per_target=1, NUM_BATCHES = num_seq_per_target //
               # batch_size would become 0 (no score computed) if batch_size > 1.
               "--batch_size", "1"]
        if chains_to_score.strip():
            cmd += ["--pdb_path_chains", chains_to_score.strip()]
        if path_to_fasta:
            cmd += ["--path_to_fasta", path_to_fasta]
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
        if proc.returncode != 0:
            return {"success": False,
                    "error": f"ProteinMPNN failed: {proc.stderr[-800:]}"}
        stem = Path(resolved).stem
        score_dir = Path(out_dir) / "score_only"
        scored = _collect_scores(score_dir, stem, fasta_seq_names)
        if not scored:
            return {"success": False,
                    "error": f"No score files found in {score_dir}. {proc.stderr[-800:]}"}
        return {
            "success": True,
            "pdb_file": resolved,
            "scored": scored,
            "note": "score = negative log-likelihood; lower means more compatible.",
        }
    except subprocess.TimeoutExpired:
        return {"success": False, "error": "ProteinMPNN timed out (300s)."}
    except Exception as exc:  # never raise out of a public tool function
        return {"success": False, "error": f"ProteinMPNN score error: {exc}"}
    finally:
        shutil.rmtree(out_dir, ignore_errors=True)


if __name__ == "__main__":
    r = proteinmpnn_design_sequences("8EK5.pdb",
            chains_to_design="E", num_sequences=2)
    print(json.dumps(r, indent=2)[:1500])
    # also exercise score-only on the same backbone
    if r.get("success"):
        s = proteinmpnn_score_sequence("8EK5.pdb",
                sequence=(r["sequences"][0]["sequence"] if r["sequences"] else ""),
                chains_to_score="E")
        print(json.dumps(s, indent=2)[:1000])
