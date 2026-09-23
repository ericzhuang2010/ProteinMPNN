# Milestone 4 - Local ProteinMPNN design for two 10LH loops

## What was requested

The September 22 email asks for Milestone 4 from the ProteinMPNN assignment and
provides two 10LH loop sequences to run through ProteinMPNN:

- `WDSSRG`
- `YTYFLD`

Milestone 4 offers a fixed-position feature. For this use case, the intended
operation is local design: the two loop segments are allowed to change and all
other residues are fixed. This follows the repository's
`examples/submit_example_4_non_fixed.sh` workflow (`--specify_non_fixed`). The
assignment PDF shortens the fixed-position example's name to `example_4.sh`;
the corresponding file in this checkout is `examples/submit_example_4.sh`.

The supplied `8EK5.pdb` is the engineered 10LH complex. The 10LH scFv is chain
E. The sequence mapping is:

| Loop sequence | PDB residue numbers | ProteinMPNN chain-E positions |
| --- | --- | --- |
| `WDSSRG` | E 200-205 | 88-93 |
| `YTYFLD` | E 339-344 | 227-232 |

ProteinMPNN positions are 1-based positions in the parsed chain, not PDB
residue numbers. The second loop is at 227-232 because the parsed chain contains
20 `X` placeholders for unresolved PDB residues E 221-240. Counting only
observed residues would incorrectly produce 207-212.

## What was implemented

`proteinmpnn_tools.py` carries forward the Milestones 2-3 design and score
functions and adds two optional, mutually exclusive arguments to
`proteinmpnn_design_sequences`:

```python
fixed_positions={"E": [1, 2, 3]}
design_only_positions={"E": [88, 89, 90]}
```

- `fixed_positions` protects the listed residues.
- `design_only_positions` allows only the listed residues to change and fixes
  every other position.

When either option is used, the wrapper invokes ProteinMPNN's official helper
scripts to parse the PDB, select designed/fixed chains, and create the
`fixed_positions_jsonl` file. It then runs `protein_mpnn_run.py` with that file.
Inputs are validated, all subprocess calls have timeouts, temporary working
files are cleaned up, and public functions retain the assignment's
JSON-serializable success/error contract.

For the requested example, chain E is designed with only these positions open:

```text
88 89 90 91 92 93 227 228 229 230 231 232
```

Chains A, B, and C are fixed context, and all other chain-E positions are
fixed.

## Delivery contents

- `proteinmpnn_tools.py` - Milestone 4 wrapper plus Milestones 2-3 behavior.
- `aibinder_tools.py` - minimal local copy of the shared PDB resolver, included
  so this delivery is independently runnable.
- `run_8ek5_two_loops.py` - reproducible 8EK5 two-loop run and validation.
- `test_proteinmpnn_tools.py` - parser, score-parser, validation, and Milestone 4
  constraint tests.
- `mpnn_out/seqs/8EK5.fa` - generated native record and eight candidate designs.
- `example_result.json` - structured result returned by the wrapper.
- `validation_report.json` - per-sample proof that no changes occurred outside
  the requested loop positions.

## Reproduce the run

From the ProteinMPNN repository root:

```bash
python assignments/06142026/deliveries/milestone4/run_8ek5_two_loops.py
```

Run the unit tests with:

```bash
python -m pytest -q \
  assignments/06142026/deliveries/milestone4/test_proteinmpnn_tools.py
```

## Verification status

The delivery was tested and run against the supplied structure:

- Unit tests: **16 passed**.
- ProteinMPNN run: **successful**, using model `v_48_020`, seed 37,
  temperature 0.1, and eight generated candidates.
- Native constrained-design score: **1.3742**.
- Lowest generated score: **0.6462** (sample 2). ProteinMPNN scores are
  negative log-likelihoods, so lower is better according to the model.
- Mutation-boundary validation: **8/8 passed**. Every sequence difference was
  inside positions 88-93 or 227-232; there were no changes outside the two
  requested loops.

The lowest-scoring generated loop pair was:

| Region | Native | Sample 2 |
| --- | --- | --- |
| Chain E positions 88-93 | `WDSSRG` | `WDTTKG` |
| Chain E positions 227-232 | `YTYFLD` | `YSWAID` |

These scores rank candidates under ProteinMPNN's sequence-backbone model; they
do not by themselves establish binding, folding, expression, or experimental
activity.
