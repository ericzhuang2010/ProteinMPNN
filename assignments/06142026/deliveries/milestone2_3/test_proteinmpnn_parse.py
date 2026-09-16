"""Unit tests for the pure FASTA parser in proteinmpnn_tools.

These tests exercise _parse_mpnn_fasta on a real-shaped ProteinMPNN output
WITHOUT importing torch or running the model, so they are fast and hermetic.
Run with: pytest tests/test_proteinmpnn_parse.py
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from proteinmpnn_tools import _parse_mpnn_fasta  # noqa: E402


# A representative ProteinMPNN output FASTA. The native header deliberately
# includes bracketed chain lists (commas INSIDE the value) to prove the parser
# does not split those on the wrong commas.
SAMPLE_FA = """\
>8EK5, score=1.0241, global_score=1.1180, fixed_chains=['A', 'B', 'C'], designed_chains=['E'], model_name=v_48_020, seed=37
EVQLVESGGGLVQPGGSLRLSCAAS
>T=0.1, sample=1, score=0.8123, global_score=0.9450, seq_recovery=0.6010
EVQLQESGPGLVKPSETLSLTCTVS
>T=0.1, sample=2, score=0.8367, global_score=0.9588, seq_recovery=0.5874
EVQLVQSGAEVKKPGASVKVSCKAS
"""


def test_record_count():
    records = _parse_mpnn_fasta(SAMPLE_FA)
    assert len(records) == 3  # native + 2 designs


def test_native_record():
    native = _parse_mpnn_fasta(SAMPLE_FA)[0]
    assert native["sequence"] == "EVQLVESGGGLVQPGGSLRLSCAAS"
    h = native["header"]
    assert h["name"] == "8EK5"
    assert h["score"] == 1.0241
    assert h["global_score"] == 1.1180
    assert h["model_name"] == "v_48_020"
    assert h["seed"] == 37
    # Bracketed chain lists must be parsed as real lists, not split on commas.
    assert h["fixed_chains"] == ["A", "B", "C"]
    assert h["designed_chains"] == ["E"]


def test_native_score_is_float():
    h = _parse_mpnn_fasta(SAMPLE_FA)[0]["header"]
    assert isinstance(h["score"], float)
    assert isinstance(h["global_score"], float)


def test_design_records():
    records = _parse_mpnn_fasta(SAMPLE_FA)
    d1, d2 = records[1], records[2]

    assert d1["sequence"] == "EVQLQESGPGLVKPSETLSLTCTVS"
    assert d1["header"]["sample"] == 1
    assert d1["header"]["T"] == 0.1
    assert d1["header"]["score"] == 0.8123
    assert d1["header"]["global_score"] == 0.9450
    assert d1["header"]["seq_recovery"] == 0.6010
    # A design record carries no structure "name" token.
    assert "name" not in d1["header"]

    assert d2["sequence"] == "EVQLVQSGAEVKKPGASVKVSCKAS"
    assert d2["header"]["sample"] == 2
    assert d2["header"]["score"] == 0.8367


def test_numeric_types_are_native():
    for rec in _parse_mpnn_fasta(SAMPLE_FA)[1:]:
        assert isinstance(rec["header"]["sample"], int)
        assert isinstance(rec["header"]["score"], float)
        assert isinstance(rec["header"]["seq_recovery"], float)


def test_multiline_sequence_is_joined():
    fa = (
        ">T=0.1, sample=1, score=0.5, global_score=0.6, seq_recovery=0.5\n"
        "EVQL\nQESG\nPGLV\n"
    )
    rec = _parse_mpnn_fasta(fa)[0]
    assert rec["sequence"] == "EVQLQESGPGLV"


def test_empty_input():
    assert _parse_mpnn_fasta("") == []
    assert _parse_mpnn_fasta("   \n  ") == []


if __name__ == "__main__":
    import pytest

    raise SystemExit(pytest.main([__file__, "-v"]))
