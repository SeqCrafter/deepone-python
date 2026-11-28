import itertools
import json
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pandas as pd
import pytest

from DeepOneCalculator import SUPPORTED_CAS9_TYPES, calculate_deepone
from DeepOneCalculator.features import _encode, _energy, _find_guides
from DeepOneCalculator.inference import _predict
from DeepOneCalculator.resources import _normalize_cas9_types

SPACER = "AT" * 10
NGG = "AAAA" + SPACER + "AGG" + "AAA"
NGN = "AAAA" + SPACER + "AGA" + "AAA"


def test_variant_names() -> None:
    assert _normalize_cas9_types("all") == list(SUPPORTED_CAS9_TYPES)
    assert _normalize_cas9_types(["SpCas9_NG", "spg", "spcas9-ng"]) == ["SpCas9-NG", "SpG"]


@pytest.mark.parametrize("value", [[], "SpRY", ["all", "SpG"], [42], None])
def test_invalid_variant_selection(value: object) -> None:
    with pytest.raises(ValueError, match="cas9_types|Cas9 type"):
        calculate_deepone(NGG, cas9_types=value)  # type: ignore[arg-type]


@pytest.mark.parametrize("pam", ["".join(b) for b in itertools.product("ATCG", repeat=3)])
def test_pam_rules_on_both_strands(pam: str) -> None:
    context = "AAAA" + SPACER + pam + "AAA"
    reverse = context.translate(str.maketrans("ATCG", "TAGC"))[::-1]
    for sequence, strand in ((context, "p"), (reverse, "m")):
        for allow_ngn in (False, True):
            guides = _find_guides(sequence, "sequence_1", allow_ngn=allow_ngn)
            expected = pam[1] == "G" and (allow_ngn or pam[2] == "G")
            assert len(guides) == int(expected)
            if expected:
                assert guides[0].context == context
                assert guides[0].grna_id == f"sequence_1_{strand}1"


def test_input_order_duplicates_and_variant_model_scope() -> None:
    sequences = [NGG] * 12  # IDs 10–12 must follow 9, not be sorted as strings.
    with patch("DeepOneCalculator.api._predict",
               side_effect=lambda x, *args: np.arange(len(x))) as predict, \
         patch("DeepOneCalculator.api._resolve_device", return_value="/CPU:0"):
        result = calculate_deepone(sequences, ["K562", "HEK293"],
                             cas9_types=["SpG", "SpCas9", "SpCas9_NG", "SpG"])
    expected_ids = [f"sequence_{index}" for index in range(1, 13)]
    assert result.sequence_id.tolist() == expected_ids * 4
    assert result.gRNA_id.nunique() == 12
    assert result.cas9_type.tolist() == ["SpG"] * 12 + ["SpCas9"] * 24 + ["SpCas9-NG"] * 12
    assert result.loc[result.cas9_type == "SpCas9", "cell_line"].tolist() == (
        ["K562"] * 12 + ["HEK293"] * 12
    )
    assert result.loc[result.cas9_type != "SpCas9", "cell_line"].isna().all()
    assert [call.args[2] for call in predict.call_args_list] == [
        "SpG", "K562", "HEK293", "SpCas9-NG",
    ]


def test_ngn_only_rows_and_missing_input_ids() -> None:
    with patch("DeepOneCalculator.api._predict",
               side_effect=lambda x, *args: np.zeros(len(x))) as predict, \
         patch("DeepOneCalculator.api._resolve_device", return_value="/CPU:0"):
        result = calculate_deepone(["A" * 30, NGN], cas9_types="all")
    assert result.sequence_id.tolist() == ["sequence_2", "sequence_2"]
    assert result.cas9_type.tolist() == ["SpCas9-NG", "SpG"]
    assert result.cell_line.isna().all()
    assert [call.args[2] for call in predict.call_args_list] == ["SpCas9-NG", "SpG"]


@pytest.mark.integration
def test_variant_reference_features_and_predictions() -> None:
    reference = json.loads(Path(__file__).with_name("variant_reference.json").read_text())
    contexts = reference["contexts"]
    energies = np.array([[_energy(s[4:27])] for s in contexts], dtype=np.float32)
    np.testing.assert_allclose(energies[:, 0], reference["energies"], rtol=1e-6, atol=1e-6)
    for model_key in ("SpCas9-NG", "SpG"):
        actual = _predict(_encode(contexts), energies, model_key, "/CPU:0", 3)
        np.testing.assert_allclose(actual, reference["scores"][model_key], rtol=1e-5, atol=1e-5)


@pytest.mark.integration
def test_real_variant_api_and_batching() -> None:
    result = calculate_deepone([NGG, NGN], ["HEK293", "K562"], cas9_types="all", device="cpu",
                         batch_size=1)
    assert len(result) == 6
    assert np.isfinite(result.deepone_score).all()
    assert result.cas9_type.tolist() == ["SpCas9", "SpCas9", "SpCas9-NG", "SpCas9-NG",
                                        "SpG", "SpG"]
    variants = calculate_deepone([NGG, NGN], "iPSC", cas9_types=["SpCas9-NG", "SpG"],
                           device="cpu", batch_size=256)
    pd.testing.assert_frame_equal(
        result.loc[result.cas9_type != "SpCas9"].reset_index(drop=True), variants,
        check_exact=False, rtol=1e-5, atol=1e-5,
    )
