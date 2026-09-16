# npz score-reading variant (reference copy)

These are **reference copies**, not part of the active codebase or test suite.

They capture the variant of `proteinmpnn_score_sequence` that reads the
score-only result from ProteinMPNN's `score_only/*.npz` files (full precision)
instead of parsing stdout.

- `proteinmpnn_tools.py` — full module with the npz-reading helpers
  (`_read_fasta_names`, `_score_from_npz`, `_collect_scores`) and a `numpy`
  dependency in `_check_deps()`.
- `npz_score_reference_tests.py` — the unit tests for the npz helpers (renamed
  so pytest does not auto-collect it; it imports symbols that only exist in this
  variant, not in the active `proteinmpnn_tools.py`).

The active code was reverted to the **stdout** approach because the assignment
asks for the "simpler reliable" option, and the npz route is more code while
being only marginally more robust (its real advantage is precision, not
reliability). See `docs/score_output_choice.md` for the full comparison and
decision.
