"""Agent-callable wrappers for ProteinMPNN design and scoring.

Milestone 4 adds residue-level constraints to the Milestones 2-3 wrapper.  A
caller can either fix listed positions or, more conveniently for local loop
design, list the only positions that ProteinMPNN may change.
"""

import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Tuple


PROTEINMPNN_DIR = os.environ.get("PROTEINMPNN_DIR", "/opt/ProteinMPNN")
RUN_SCRIPT = os.path.join(PROTEINMPNN_DIR, "protein_mpnn_run.py")
HELPER_DIR = os.path.join(PROTEINMPNN_DIR, "helper_scripts")
PARSE_CHAINS_SCRIPT = os.path.join(HELPER_DIR, "parse_multiple_chains.py")
ASSIGN_CHAINS_SCRIPT = os.path.join(HELPER_DIR, "assign_fixed_chains.py")
MAKE_FIXED_SCRIPT = os.path.join(HELPER_DIR, "make_fixed_positions_dict.py")

_TORCH_OK = False
try:
    import torch  # noqa: F401

    _TORCH_OK = True
except ImportError:
    pass

from aibinder_tools import _find_pdb  # noqa: E402


def _check_deps(require_position_helpers: bool = False) -> Optional[Dict[str, Any]]:
    """Return an error dictionary when a required dependency is unavailable."""
    missing = []
    if not _TORCH_OK:
        missing.append("torch")
    if not os.path.isfile(RUN_SCRIPT):
        missing.append(
            f"ProteinMPNN (set PROTEINMPNN_DIR; looked in {PROTEINMPNN_DIR})"
        )
    if require_position_helpers:
        for label, path in (
            ("parse_multiple_chains.py", PARSE_CHAINS_SCRIPT),
            ("assign_fixed_chains.py", ASSIGN_CHAINS_SCRIPT),
            ("make_fixed_positions_dict.py", MAKE_FIXED_SCRIPT),
        ):
            if not os.path.isfile(path):
                missing.append(f"ProteinMPNN helper {label}")
    if missing:
        return {
            "success": False,
            "error": (
                "Missing dependencies: "
                + ", ".join(missing)
                + ". Install torch and set PROTEINMPNN_DIR to the ProteinMPNN checkout."
            ),
        }
    return None


_FLOAT_KEYS = {"score", "global_score", "seq_recovery", "T"}
_INT_KEYS = {"sample", "seed"}
_FIELD_SPLIT_RE = re.compile(r",(?![^\[]*\])")


def _parse_header_value(key: str, value: str) -> Any:
    """Convert a ProteinMPNN FASTA-header value to a native Python value."""
    value = value.strip()
    if value.startswith("["):
        return [
            part.strip().strip("'\"")
            for part in value.strip("[]").split(",")
            if part.strip()
        ]
    if key in _FLOAT_KEYS:
        try:
            return float(value)
        except ValueError:
            return value
    if key in _INT_KEYS:
        try:
            return int(value)
        except ValueError:
            return value
    return value


def _parse_mpnn_fasta(text: str) -> List[Dict[str, Any]]:
    """Parse ProteinMPNN output FASTA into JSON-serializable records."""
    records = []
    blocks = [block for block in text.strip().split(">") if block.strip()]
    for block in blocks:
        lines = block.splitlines()
        header_line = lines[0]
        sequence = "".join(line.strip() for line in lines[1:])
        metadata: Dict[str, Any] = {}
        for index, piece in enumerate(_FIELD_SPLIT_RE.split(header_line)):
            piece = piece.strip()
            if "=" in piece:
                key, value = piece.split("=", 1)
                key = key.strip()
                metadata[key] = _parse_header_value(key, value)
            elif index == 0 and piece:
                metadata["name"] = piece
        records.append({"header": metadata, "sequence": sequence})
    return records


def _normalize_position_mapping(
    positions: Optional[Mapping[str, List[int]]], label: str
) -> Dict[str, List[int]]:
    """Validate and normalize a chain-to-1-based-position mapping.

    Internal helpers may raise ``ValueError``; public functions catch these and
    return the standard ``{"success": False, "error": ...}`` shape.
    """
    if positions is None:
        return {}
    if not isinstance(positions, Mapping):
        raise ValueError(f"{label} must be a mapping such as {{'E': [88, 89]}}.")

    normalized: Dict[str, List[int]] = {}
    for raw_chain, raw_positions in positions.items():
        chain = str(raw_chain).strip()
        if not chain or any(character.isspace() for character in chain):
            raise ValueError(f"Invalid chain identifier in {label}: {raw_chain!r}.")
        if isinstance(raw_positions, (str, bytes)) or not isinstance(
            raw_positions, (list, tuple, set)
        ):
            raise ValueError(f"{label}[{chain!r}] must be a list of positions.")

        cleaned = []
        for raw_position in raw_positions:
            if isinstance(raw_position, bool):
                raise ValueError(f"Positions in {label} must be positive integers.")
            try:
                position = int(raw_position)
            except (TypeError, ValueError) as exc:
                raise ValueError(
                    f"Positions in {label} must be positive integers."
                ) from exc
            if position < 1 or position != raw_position:
                raise ValueError(f"Positions in {label} must be positive integers.")
            cleaned.append(position)
        if not cleaned:
            raise ValueError(f"{label}[{chain!r}] cannot be empty.")
        normalized[chain] = sorted(set(cleaned))

    if not normalized:
        raise ValueError(f"{label} cannot be empty when provided.")
    return normalized


def _resolve_designed_chains(
    chains_to_design: str, position_mapping: Mapping[str, List[int]]
) -> List[str]:
    """Resolve chain order and ensure every constrained design chain is explicit."""
    chains = [chain for chain in chains_to_design.split() if chain]
    if position_mapping:
        if not chains:
            chains = list(position_mapping)
        if set(chains) != set(position_mapping):
            raise ValueError(
                "When residue constraints are used, their chain keys must exactly "
                "match chains_to_design."
            )
    return chains


def _position_list_argument(
    chains: List[str], positions: Mapping[str, List[int]]
) -> str:
    """Format positions for make_fixed_positions_dict.py in chain-list order."""
    return ", ".join(" ".join(str(value) for value in positions[chain]) for chain in chains)


def _run_checked(command: List[str], timeout: int = 300) -> subprocess.CompletedProcess:
    """Run one ProteinMPNN/helper command or raise a concise runtime error."""
    process = subprocess.run(
        command, capture_output=True, text=True, timeout=timeout, check=False
    )
    if process.returncode != 0:
        detail = (process.stderr or process.stdout or "no diagnostic output")[-1200:]
        raise RuntimeError(f"command failed ({Path(command[1]).name}): {detail}")
    return process


def _prepare_position_constraints(
    resolved_pdb: str,
    workspace: Path,
    designed_chains: List[str],
    positions: Mapping[str, List[int]],
    design_only: bool,
) -> Tuple[Path, Path, Path]:
    """Build parsed-chain, assigned-chain, and fixed-position JSONL files."""
    pdb_dir = workspace / "pdbs"
    pdb_dir.mkdir()
    shutil.copy2(resolved_pdb, pdb_dir / Path(resolved_pdb).name)

    parsed_path = workspace / "parsed_pdbs.jsonl"
    assigned_path = workspace / "assigned_pdbs.jsonl"
    fixed_path = workspace / "fixed_pdbs.jsonl"
    chain_argument = " ".join(designed_chains)

    _run_checked(
        [
            sys.executable,
            PARSE_CHAINS_SCRIPT,
            "--input_path",
            str(pdb_dir),
            "--output_path",
            str(parsed_path),
        ]
    )
    _run_checked(
        [
            sys.executable,
            ASSIGN_CHAINS_SCRIPT,
            "--input_path",
            str(parsed_path),
            "--output_path",
            str(assigned_path),
            "--chain_list",
            chain_argument,
        ]
    )
    fixed_command = [
        sys.executable,
        MAKE_FIXED_SCRIPT,
        "--input_path",
        str(parsed_path),
        "--output_path",
        str(fixed_path),
        "--chain_list",
        chain_argument,
        "--position_list",
        _position_list_argument(designed_chains, positions),
    ]
    if design_only:
        fixed_command.append("--specify_non_fixed")
    _run_checked(fixed_command)
    return parsed_path, assigned_path, fixed_path


def _changed_positions(native_sequence: str, designed_sequence: str) -> List[int]:
    """Return 1-based positions that differ between equal-length sequences."""
    if len(native_sequence) != len(designed_sequence):
        raise ValueError("Native and designed sequences must have equal length.")
    return [
        index
        for index, (native, designed) in enumerate(
            zip(native_sequence, designed_sequence), start=1
        )
        if native != designed
    ]


def proteinmpnn_design_sequences(
    pdb_path: str,
    chains_to_design: str = "",
    num_sequences: int = 8,
    sampling_temp: float = 0.1,
    seed: int = 37,
    model_name: str = "v_48_020",
    fixed_positions: Optional[Mapping[str, List[int]]] = None,
    design_only_positions: Optional[Mapping[str, List[int]]] = None,
) -> Dict[str, Any]:
    """Design sequences for a backbone, optionally with residue constraints.

    ``fixed_positions`` lists residues that must retain their native identity.
    ``design_only_positions`` lists the only residues that may change and fixes
    the complement. Positions are 1-based indices in ProteinMPNN's parsed chain,
    not PDB residue numbers. The two options are mutually exclusive.

    Score is negative log-likelihood; lower means more compatible. The public
    function always returns a success/error dictionary and never raises.
    """
    out_dir: Optional[str] = None
    try:
        if fixed_positions is not None and design_only_positions is not None:
            return {
                "success": False,
                "error": "Use fixed_positions or design_only_positions, not both.",
            }

        position_label = (
            "design_only_positions"
            if design_only_positions is not None
            else "fixed_positions"
        )
        raw_positions = (
            design_only_positions
            if design_only_positions is not None
            else fixed_positions
        )
        normalized_positions = _normalize_position_mapping(raw_positions, position_label)
        designed_chains = _resolve_designed_chains(
            chains_to_design, normalized_positions
        )

        dependency_error = _check_deps(
            require_position_helpers=bool(normalized_positions)
        )
        if dependency_error:
            return dependency_error

        resolved = _find_pdb(pdb_path)
        if not resolved:
            return {"success": False, "error": f"PDB file not found: {pdb_path}"}

        out_dir = tempfile.mkdtemp(prefix="mpnn_")
        workspace = Path(out_dir)
        model_output = workspace / "model_output"

        command = [
            sys.executable,
            RUN_SCRIPT,
            "--out_folder",
            str(model_output),
            "--num_seq_per_target",
            str(num_sequences),
            "--sampling_temp",
            str(sampling_temp),
            "--seed",
            str(seed),
            "--model_name",
            model_name,
            "--batch_size",
            "1",
        ]

        if normalized_positions:
            parsed, assigned, fixed = _prepare_position_constraints(
                resolved,
                workspace,
                designed_chains,
                normalized_positions,
                design_only=design_only_positions is not None,
            )
            command.extend(
                [
                    "--jsonl_path",
                    str(parsed),
                    "--chain_id_jsonl",
                    str(assigned),
                    "--fixed_positions_jsonl",
                    str(fixed),
                ]
            )
        else:
            command.extend(["--pdb_path", resolved])
            if designed_chains:
                command.extend(["--pdb_path_chains", " ".join(designed_chains)])

        process = _run_checked(command)
        fasta_path = model_output / "seqs" / f"{Path(resolved).stem}.fa"
        if not fasta_path.is_file():
            return {
                "success": False,
                "error": (
                    f"ProteinMPNN produced no FASTA at {fasta_path}. "
                    f"{process.stderr[-800:]}"
                ),
            }

        records = _parse_mpnn_fasta(fasta_path.read_text())
        if not records:
            return {"success": False, "error": f"Could not parse {fasta_path}."}

        kept_fasta = tempfile.NamedTemporaryFile(
            prefix=f"{Path(resolved).stem}_mpnn_", suffix=".fa", delete=False
        )
        kept_fasta.close()
        shutil.copyfile(fasta_path, kept_fasta.name)

        native = records[0]
        native_header = native["header"]
        sequences = []
        for design in records[1:]:
            header = design["header"]
            record = {
                "sample": header.get("sample"),
                "sequence": design["sequence"],
                "score": header.get("score"),
                "global_score": header.get("global_score"),
                "seq_recovery": header.get("seq_recovery"),
                "temperature": header.get("T"),
            }
            if design_only_positions is not None and len(designed_chains) == 1:
                record["changed_positions"] = _changed_positions(
                    native["sequence"], design["sequence"]
                )
            sequences.append(record)

        result: Dict[str, Any] = {
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
            "out_fasta_path": kept_fasta.name,
        }
        if normalized_positions:
            result["position_constraint_mode"] = (
                "design_only" if design_only_positions is not None else "fixed"
            )
            result[position_label] = normalized_positions
        return result
    except subprocess.TimeoutExpired:
        return {"success": False, "error": "ProteinMPNN timed out (300s)."}
    except Exception as exc:
        return {"success": False, "error": f"ProteinMPNN design error: {exc}"}
    finally:
        if out_dir:
            shutil.rmtree(out_dir, ignore_errors=True)


_SCORE_RE = re.compile(
    r"Score for (?P<name>.+?) from (?P<src>PDB|FASTA), mean: "
    r"(?P<score>[-\d.eE]+),.*?global score, mean: (?P<gscore>[-\d.eE]+)"
)


def _parse_score_stdout(text: str) -> List[Dict[str, Any]]:
    """Parse ProteinMPNN ``--score_only`` stdout."""
    return [
        {
            "name": match.group("name").strip(),
            "score": float(match.group("score")),
            "global_score": float(match.group("gscore")),
        }
        for match in _SCORE_RE.finditer(text)
    ]


def proteinmpnn_score_sequence(
    pdb_path: str,
    sequence: str = "",
    fasta_path: str = "",
    chains_to_score: str = "",
    seed: int = 37,
) -> Dict[str, Any]:
    """Score sequence(s) against a PDB backbone; lower scores are better."""
    dependency_error = _check_deps()
    if dependency_error:
        return dependency_error
    resolved = _find_pdb(pdb_path)
    if not resolved:
        return {"success": False, "error": f"PDB file not found: {pdb_path}"}

    out_dir = tempfile.mkdtemp(prefix="mpnn_")
    try:
        input_fasta = ""
        if fasta_path:
            if not os.path.isfile(fasta_path):
                return {
                    "success": False,
                    "error": f"FASTA file not found: {fasta_path}",
                }
            input_fasta = os.path.abspath(fasta_path)
        elif sequence.strip():
            temporary_fasta = Path(out_dir) / "input_to_score.fa"
            temporary_fasta.write_text(f">seq1\n{sequence.strip()}\n")
            input_fasta = str(temporary_fasta)

        command = [
            sys.executable,
            RUN_SCRIPT,
            "--pdb_path",
            resolved,
            "--out_folder",
            out_dir,
            "--score_only",
            "1",
            "--seed",
            str(seed),
            "--num_seq_per_target",
            "1",
            "--batch_size",
            "1",
        ]
        if chains_to_score.strip():
            command.extend(["--pdb_path_chains", chains_to_score.strip()])
        if input_fasta:
            command.extend(["--path_to_fasta", input_fasta])

        process = _run_checked(command)
        scored = _parse_score_stdout(process.stdout)
        if not scored:
            return {
                "success": False,
                "error": f"No scores parsed from ProteinMPNN output. {process.stdout[-800:]}",
            }
        return {
            "success": True,
            "pdb_file": resolved,
            "scored": scored,
            "note": "score = negative log-likelihood; lower means more compatible.",
        }
    except subprocess.TimeoutExpired:
        return {"success": False, "error": "ProteinMPNN timed out (300s)."}
    except Exception as exc:
        return {"success": False, "error": f"ProteinMPNN score error: {exc}"}
    finally:
        shutil.rmtree(out_dir, ignore_errors=True)


if __name__ == "__main__":
    result = proteinmpnn_design_sequences(
        "8EK5.pdb",
        chains_to_design="E",
        num_sequences=2,
        design_only_positions={"E": [88, 89, 90, 91, 92, 93, 227, 228, 229, 230, 231, 232]},
    )
    print(json.dumps(result, indent=2)[:3000])
