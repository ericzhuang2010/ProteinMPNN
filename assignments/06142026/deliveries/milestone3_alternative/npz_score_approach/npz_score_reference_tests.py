"""Unit tests for the score-only .npz reading in proteinmpnn_tools.

These tests build fake ProteinMPNN score_only/*.npz files with numpy and feed
them to _collect_scores / _score_from_npz, so they exercise the score-parsing
logic WITHOUT running the model. Run with: pytest tests/test_proteinmpnn_score.py
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np  # noqa: E402

from proteinmpnn_tools import _collect_scores, _score_from_npz  # noqa: E402


def _write_npz(path, score, global_score):
    """Write a score_only-style .npz with `score` and `global_score` arrays."""
    np.savez(str(path), score=np.array(score, dtype=np.float32),
             global_score=np.array(global_score, dtype=np.float32))


def test_score_from_npz_averages_and_is_native(tmp_path):
    f = tmp_path / "x_pdb.npz"
    _write_npz(f, score=[1.0, 2.0], global_score=[3.0, 5.0])
    rec = _score_from_npz(f, "mylabel")
    assert rec["name"] == "mylabel"
    # arrays are averaged: mean([1.0, 2.0]) == 1.5, mean([3.0, 5.0]) == 4.0
    assert rec["score"] == 1.5
    assert rec["global_score"] == 4.0
    # ...and cast to native python float (JSON-serializable, not numpy).
    assert type(rec["score"]) is float
    assert type(rec["global_score"]) is float


def test_collect_scores_pdb_only(tmp_path):
    score_dir = tmp_path / "score_only"
    score_dir.mkdir()
    _write_npz(score_dir / "8EK5_pdb.npz", score=[1.1427], global_score=[1.4197])

    scored = _collect_scores(score_dir, "8EK5", fasta_seq_names=[])
    assert len(scored) == 1
    assert scored[0]["name"] == "8EK5 (native from PDB)"
    assert round(scored[0]["score"], 4) == 1.1427
    assert round(scored[0]["global_score"], 4) == 1.4197


def test_collect_scores_pdb_plus_fasta_with_names(tmp_path):
    score_dir = tmp_path / "score_only"
    score_dir.mkdir()
    _write_npz(score_dir / "8EK5_pdb.npz", score=[1.1427], global_score=[1.4197])
    _write_npz(score_dir / "8EK5_fasta_1.npz", score=[0.7068], global_score=[1.26])
    _write_npz(score_dir / "8EK5_fasta_2.npz", score=[0.6798], global_score=[1.25])

    scored = _collect_scores(score_dir, "8EK5", fasta_seq_names=["designA", "designB"])
    assert [s["name"] for s in scored] == [
        "8EK5 (native from PDB)", "designA", "designB",
    ]
    assert round(scored[1]["score"], 4) == 0.7068
    assert round(scored[2]["score"], 4) == 0.6798


def test_collect_scores_fasta_name_fallback(tmp_path):
    # More .npz files than provided names -> fall back to "fasta_<fc>".
    score_dir = tmp_path / "score_only"
    score_dir.mkdir()
    _write_npz(score_dir / "S_pdb.npz", score=[1.0], global_score=[1.0])
    _write_npz(score_dir / "S_fasta_1.npz", score=[0.9], global_score=[1.0])
    _write_npz(score_dir / "S_fasta_2.npz", score=[0.8], global_score=[1.0])

    scored = _collect_scores(score_dir, "S", fasta_seq_names=["only_one"])
    assert [s["name"] for s in scored] == [
        "S (native from PDB)", "only_one", "fasta_2",
    ]


def test_collect_scores_empty_when_no_files(tmp_path):
    score_dir = tmp_path / "score_only"
    score_dir.mkdir()
    assert _collect_scores(score_dir, "nothing", fasta_seq_names=[]) == []


if __name__ == "__main__":
    import pytest

    raise SystemExit(pytest.main([__file__, "-v"]))
