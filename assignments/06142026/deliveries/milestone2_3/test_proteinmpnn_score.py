"""Tests for the score-only stdout parser and the error/edge-case contract.

All tests here run WITHOUT torch and WITHOUT invoking ProteinMPNN: the stdout
parser is a pure function, and the error paths are reached before any
subprocess call (so they're fast and hermetic). They cover the Definition of
Done items: missing-dependency & bad-input cases return an error dict and never
crash, and every returned value is JSON-serializable.

Run with: pytest tests/test_proteinmpnn_score.py
"""

import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import proteinmpnn_tools as t  # noqa: E402
from proteinmpnn_tools import _parse_score_stdout  # noqa: E402


# A realistic chunk of ProteinMPNN --score_only stdout: setup chatter plus the
# two summary lines we care about (one PDB native, one external FASTA sequence).
SAMPLE_STDOUT = """\
----------------------------------------
Number of edges: 48
Training noise level: 0.2A
chain_id_jsonl is NOT loaded
Score for 8EK5 from PDB, mean: 1.1427, std: 0.0000, sample size: 1,  global score, mean: 1.4197, std: 0.0000, sample size: 1
Score for 8EK5_1 from FASTA, mean: 0.7068, std: 0.0000, sample size: 1,  global score, mean: 1.2600, std: 0.0000, sample size: 1
"""


def test_parse_score_stdout_parses_pdb_and_fasta():
    scored = _parse_score_stdout(SAMPLE_STDOUT)
    assert len(scored) == 2
    assert scored[0] == {"name": "8EK5", "score": 1.1427, "global_score": 1.4197}
    assert scored[1] == {"name": "8EK5_1", "score": 0.7068, "global_score": 1.2600}


def test_parse_score_stdout_values_are_native_floats():
    for rec in _parse_score_stdout(SAMPLE_STDOUT):
        assert type(rec["score"]) is float
        assert type(rec["global_score"]) is float


def test_parse_score_stdout_ignores_non_score_lines():
    noise = "loading model...\nseed=37\nDone in 1.2s\n"
    assert _parse_score_stdout(noise) == []


def test_parse_score_stdout_empty():
    assert _parse_score_stdout("") == []


def test_missing_torch_returns_error_dict(monkeypatch):
    # Simulate torch being unavailable -> _check_deps reports it, never crashes.
    monkeypatch.setattr(t, "_TORCH_OK", False)
    result = t.proteinmpnn_design_sequences("8EK5.pdb", chains_to_design="E")
    assert result["success"] is False
    assert "torch" in result["error"]
    json.dumps(result)  # must be JSON-serializable


def test_missing_proteinmpnn_returns_error_dict(monkeypatch):
    # torch present, but the run script can't be found.
    monkeypatch.setattr(t, "_TORCH_OK", True)
    monkeypatch.setattr(t, "RUN_SCRIPT", "/nonexistent/protein_mpnn_run.py")
    result = t.proteinmpnn_score_sequence("8EK5.pdb", sequence="AAAA")
    assert result["success"] is False
    assert "ProteinMPNN" in result["error"]
    json.dumps(result)


def test_bad_pdb_returns_error_dict_design(monkeypatch):
    # Deps satisfied (point RUN_SCRIPT at any real file) so we reach the PDB
    # lookup, then pass a path that cannot resolve -> clean error, no subprocess.
    monkeypatch.setattr(t, "_TORCH_OK", True)
    monkeypatch.setattr(t, "RUN_SCRIPT", __file__)
    result = t.proteinmpnn_design_sequences("definitely_missing_12345.pdb")
    assert result["success"] is False
    assert "PDB file not found" in result["error"]
    json.dumps(result)


def test_bad_pdb_returns_error_dict_score(monkeypatch):
    monkeypatch.setattr(t, "_TORCH_OK", True)
    monkeypatch.setattr(t, "RUN_SCRIPT", __file__)
    result = t.proteinmpnn_score_sequence("definitely_missing_12345.pdb", sequence="AAAA")
    assert result["success"] is False
    assert "PDB file not found" in result["error"]
    json.dumps(result)


if __name__ == "__main__":
    import pytest

    raise SystemExit(pytest.main([__file__, "-v"]))
