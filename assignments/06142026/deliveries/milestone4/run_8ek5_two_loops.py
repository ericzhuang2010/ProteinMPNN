"""Reproduce the Milestone 4 local-design example on the supplied 8EK5 PDB."""

import json
import os
import shutil
from pathlib import Path


DELIVERY_DIR = Path(__file__).resolve().parent
REPO_ROOT = DELIVERY_DIR.parents[3]
PDB_PATH = REPO_ROOT / "assignments/06142026/requirements/8EK5.pdb"

# These must be set before proteinmpnn_tools is imported because it resolves
# tool paths at module-import time.
os.environ.setdefault("PROTEINMPNN_DIR", str(REPO_ROOT))
os.environ.setdefault("AIBINDER_BASE_DIR", str(PDB_PATH.parent))

from proteinmpnn_tools import proteinmpnn_design_sequences  # noqa: E402


DESIGN_ONLY_POSITIONS = {
    "E": [88, 89, 90, 91, 92, 93, 227, 228, 229, 230, 231, 232]
}


def main() -> int:
    """Run ProteinMPNN, retain the FASTA, and write machine-readable validation."""
    result = proteinmpnn_design_sequences(
        str(PDB_PATH),
        chains_to_design="E",
        num_sequences=8,
        sampling_temp=0.1,
        seed=37,
        design_only_positions=DESIGN_ONLY_POSITIONS,
    )
    if not result.get("success"):
        (DELIVERY_DIR / "example_result.json").write_text(
            json.dumps(result, indent=2) + "\n"
        )
        print(json.dumps(result, indent=2))
        return 1

    output_fasta = DELIVERY_DIR / "mpnn_out/seqs/8EK5.fa"
    output_fasta.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(result["out_fasta_path"], output_fasta)
    result["out_fasta_path"] = str(output_fasta)

    allowed = set(DESIGN_ONLY_POSITIONS["E"])
    validations = []
    for design in result["sequences"]:
        changed = design.get("changed_positions", [])
        outside = sorted(set(changed) - allowed)
        validations.append(
            {
                "sample": design["sample"],
                "changed_positions": changed,
                "outside_allowed_positions": outside,
                "valid": not outside,
            }
        )

    validation_report = {
        "success": all(item["valid"] for item in validations),
        "chain": "E",
        "design_only_positions": sorted(allowed),
        "sequence_to_position_mapping": {
            "WDSSRG": {
                "proteinmpnn_positions": [88, 89, 90, 91, 92, 93],
                "pdb_residues": "E 200-205",
            },
            "YTYFLD": {
                "proteinmpnn_positions": [227, 228, 229, 230, 231, 232],
                "pdb_residues": "E 339-344",
            },
        },
        "samples": validations,
    }
    (DELIVERY_DIR / "example_result.json").write_text(
        json.dumps(result, indent=2) + "\n"
    )
    (DELIVERY_DIR / "validation_report.json").write_text(
        json.dumps(validation_report, indent=2) + "\n"
    )
    print(json.dumps(validation_report, indent=2))
    return 0 if validation_report["success"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
