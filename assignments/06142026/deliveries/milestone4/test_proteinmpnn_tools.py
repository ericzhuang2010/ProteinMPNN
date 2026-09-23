"""Pure-function and error-contract tests for the Milestone 4 delivery."""

import json
import sys
from pathlib import Path

import pytest


DELIVERY_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(DELIVERY_DIR))

import proteinmpnn_tools as tools  # noqa: E402


SAMPLE_FASTA = """\
>8EK5, score=1.0241, global_score=1.1180, fixed_chains=['A', 'B', 'C'], designed_chains=['E'], model_name=v_48_020, seed=37
VLTQPPSVSVSPGQTA
>T=0.1, sample=1, score=0.8123, global_score=0.9450, seq_recovery=0.6010
VLTQPPSVSVSPGQAA
"""


def test_fasta_parser_retains_native_types_and_chain_lists():
    records = tools._parse_mpnn_fasta(SAMPLE_FASTA)
    assert len(records) == 2
    assert records[0]["header"]["fixed_chains"] == ["A", "B", "C"]
    assert records[0]["header"]["designed_chains"] == ["E"]
    assert records[0]["header"]["score"] == 1.0241
    assert records[1]["header"]["sample"] == 1
    assert records[1]["header"]["T"] == 0.1


def test_normalize_position_mapping_sorts_and_deduplicates():
    assert tools._normalize_position_mapping(
        {"E": [93, 88, 89, 88]}, "design_only_positions"
    ) == {"E": [88, 89, 93]}


@pytest.mark.parametrize(
    "bad",
    [
        {"E": []},
        {"E": [0]},
        {"E": [-1]},
        {"E": [1.5]},
        {"E": [True]},
        {"E": "88 89"},
    ],
)
def test_normalize_position_mapping_rejects_bad_values(bad):
    with pytest.raises(ValueError):
        tools._normalize_position_mapping(bad, "design_only_positions")


def test_constraint_chain_keys_must_match_designed_chains():
    with pytest.raises(ValueError):
        tools._resolve_designed_chains("A E", {"E": [88]})


def test_design_chain_can_be_inferred_from_constraint_mapping():
    assert tools._resolve_designed_chains("", {"E": [88]}) == ["E"]


def test_position_argument_follows_chain_order():
    positions = {"E": [88, 89], "A": [1, 2]}
    assert tools._position_list_argument(["A", "E"], positions) == "1 2, 88 89"


def test_changed_positions_are_one_based():
    assert tools._changed_positions("ABCDE", "ABXDY") == [3, 5]


def test_changed_positions_reject_different_lengths():
    with pytest.raises(ValueError):
        tools._changed_positions("ABC", "AB")


def test_mutually_exclusive_constraint_modes_return_error_dict():
    result = tools.proteinmpnn_design_sequences(
        "8EK5.pdb",
        fixed_positions={"E": [1]},
        design_only_positions={"E": [2]},
    )
    assert result["success"] is False
    assert "not both" in result["error"]
    json.dumps(result)


def test_bad_constraint_returns_error_instead_of_raising(monkeypatch):
    monkeypatch.setattr(tools, "_TORCH_OK", True)
    monkeypatch.setattr(tools, "RUN_SCRIPT", __file__)
    result = tools.proteinmpnn_design_sequences(
        "8EK5.pdb", design_only_positions={"E": [0]}
    )
    assert result["success"] is False
    assert "positive integers" in result["error"]
    json.dumps(result)


def test_score_stdout_parser():
    stdout = (
        "Score for 8EK5 from PDB, mean: 1.1427, std: 0.0000, sample size: 1, "
        " global score, mean: 1.4197, std: 0.0000, sample size: 1\n"
    )
    assert tools._parse_score_stdout(stdout) == [
        {"name": "8EK5", "score": 1.1427, "global_score": 1.4197}
    ]
