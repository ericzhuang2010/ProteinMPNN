# ProteinMPNN Milestone 1 Run Notes - 8EK5

## Scope

All Milestone 1 runs use the public RCSB structure `8EK5`. Chain `E` is the
scFv binder being designed or scored. Chains `A`, `B`, and `C` remain fixed as
the structural context.

ProteinMPNN scores are negative log-likelihood values. Lower scores indicate
greater sequence-backbone compatibility according to the model; they do not
by themselves prove folding, expression, binding, or biological function.

## Run 1: Design Four Sequences for 8EK5 Chain E

### Command

```bash
python3 protein_mpnn_run.py \
  --pdb_path 8EK5.pdb \
  --pdb_path_chains "E" \
  --out_folder ./outputs/mpnn_out \
  --num_seq_per_target 4 \
  --sampling_temp "0.1" \
  --seed 37 \
  --batch_size 1
```

The original PDB was later copied to the repository root as `8EK5.pdb`. The
FASTA result was preserved at:

`assignments/mpnn_out/design_8EK5/8EK5.fa`

### Console output

```text
chain_id_jsonl is NOT loaded
fixed_positions_jsonl is NOT loaded
pssm_jsonl is NOT loaded
omit_AA_jsonl is NOT loaded
bias_AA_jsonl is NOT loaded
tied_positions_jsonl is NOT loaded
bias by residue dictionary is not loaded, or not provided

Number of edges: 48
Training noise level: 0.2A
Generating sequences for: 8EK5
4 sequences of length 631 generated in 7.1442 seconds
```

### Saved FASTA results

| Record | Score | Global score | Sequence recovery |
| --- | ---: | ---: | ---: |
| Native 8EK5 chain E | 1.1251 | 1.4086 | N/A |
| Design 1 | 0.6835 | 1.2436 | 0.6388 |
| Design 2 | 0.6451 | 1.2276 | 0.6432 |
| Design 3 | 0.6647 | 1.2433 | 0.6344 |
| Design 4 | 0.6469 | 1.2212 | 0.6123 |

The console's length `631` is the full model input length. Each sequence in
the saved chain-E FASTA is 247 characters, including the `X` separator between
the scFv variable domains.

## Run 2: Score the Native 8EK5 Chain-E Sequence

This run uses score-only mode without an external FASTA, so ProteinMPNN scores
the native chain-E sequence present in `8EK5.pdb`.

### Command

```bash
python3 protein_mpnn_run.py \
  --pdb_path 8EK5.pdb \
  --pdb_path_chains "E" \
  --out_folder assignments/mpnn_out/score_only_8EK5_no_fasta \
  --num_seq_per_target 10 \
  --sampling_temp "0.1" \
  --score_only 1 \
  --seed 37 \
  --batch_size 1
```

### Console output

```text
Number of edges: 48
Training noise level: 0.2A
Score for 8EK5 from PDB, mean: 1.1392, std: 0.0137, sample size: 10, global score, mean: 1.4132, std: 0.0068, sample size: 10
```

### Summary

- Target: `8EK5`
- Scored chain: `E`
- Sequence source: native sequence from the PDB
- Score mean/std: `1.1392 / 0.0137`
- Global score mean/std: `1.4132 / 0.0068`
- Sample size: `10`
- Output: `assignments/mpnn_out/score_only_8EK5_no_fasta/score_only/8EK5_pdb.npz`

## Run 3: Score the Native and Four Designed 8EK5 Sequences

This run uses the FASTA produced in Run 1. The FASTA contains the native
chain-E record followed by the four designed sequences.

### Command

```bash
python3 protein_mpnn_run.py \
  --path_to_fasta assignments/mpnn_out/design_8EK5/8EK5.fa \
  --pdb_path 8EK5.pdb \
  --pdb_path_chains "E" \
  --out_folder assignments/mpnn_out/score_only_8EK5_fasta \
  --num_seq_per_target 5 \
  --sampling_temp "0.1" \
  --score_only 1 \
  --seed 13 \
  --batch_size 1
```

### Console output

```text
Number of edges: 48
Training noise level: 0.2A
Score for 8EK5 from PDB, mean: 1.1395, std: 0.0063, sample size: 5, global score, mean: 1.4156, std: 0.0060, sample size: 5
Score for 8EK5_1 from FASTA, mean: 1.1285, std: 0.0162, sample size: 5, global score, mean: 1.4121, std: 0.0082, sample size: 5
Score for 8EK5_2 from FASTA, mean: 0.7517, std: 0.0087, sample size: 5, global score, mean: 1.2685, std: 0.0094, sample size: 5
Score for 8EK5_3 from FASTA, mean: 0.6913, std: 0.0084, sample size: 5, global score, mean: 1.2431, std: 0.0060, sample size: 5
Score for 8EK5_4 from FASTA, mean: 0.7164, std: 0.0096, sample size: 5, global score, mean: 1.2557, std: 0.0084, sample size: 5
Score for 8EK5_5 from FASTA, mean: 0.6889, std: 0.0120, sample size: 5, global score, mean: 1.2465, std: 0.0062, sample size: 5
```

### Saved score-only results

| Output | FASTA mapping | Mean score | Mean global score |
| --- | --- | ---: | ---: |
| `8EK5_pdb.npz` | Native sequence read directly from PDB | 1.1395 | 1.4156 |
| `8EK5_fasta_1.npz` | Native record from FASTA | 1.1285 | 1.4121 |
| `8EK5_fasta_2.npz` | Design 1 | 0.7517 | 1.2685 |
| `8EK5_fasta_3.npz` | Design 2 | 0.6913 | **1.2431** |
| `8EK5_fasta_4.npz` | Design 3 | 0.7164 | 1.2557 |
| `8EK5_fasta_5.npz` | Design 4 | **0.6889** | 1.2465 |

The output directory is:

`assignments/mpnn_out/score_only_8EK5_fasta/score_only/`

Among the four designs, Design 4 has the lowest mean score (`0.6889`), while
Design 2 has the lowest mean global score (`1.2431`). All four designs score
better than the native 8EK5 sequence under these score-only evaluations.

## Flag Reference

| Flag | Meaning in these runs |
| --- | --- |
| `--pdb_path` | Selects the `8EK5` structure used as the backbone. |
| `--pdb_path_chains "E"` | Designs or scores chain `E`; chains `A`, `B`, and `C` provide fixed context. |
| `--out_folder` | Selects where ProteinMPNN writes FASTA or `.npz` results. |
| `--num_seq_per_target` | Requests four designs in Run 1 or sets the number of score passes in Runs 2 and 3. |
| `--sampling_temp "0.1"` | Uses a conservative sampling temperature during design. No new sequence is generated in score-only mode. |
| `--score_only 1` | Scores structure-sequence compatibility without generating new sequences. |
| `--path_to_fasta` | Supplies the native and four designed chain-E sequences for Run 3. |
| `--seed` | Makes the stochastic run reproducible. |
| `--batch_size 1` | Uses the safest CPU-friendly batch size. |

## Milestone 1 Conclusion

Milestone 1 is complete for public structure `8EK5`:

1. ProteinMPNN was run manually to design chain `E`, producing four sequences.
2. The native chain-E sequence was evaluated in score-only mode.
3. The native FASTA record and all four designs were evaluated in score-only mode.
4. The commands, flags, outputs, and score-file mappings are documented above.
