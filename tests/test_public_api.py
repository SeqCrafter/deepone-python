"""验证统一的评分入口及旧名称的移除。"""

import inspect
from importlib import import_module

import pytest

import DeepOneCalculator as calculator


@pytest.mark.parametrize(("name", "module"), [
    ("calculate_deepone", "api"),
    ("calculate_crispron", "crispron"),
    ("calculate_crisprspec", "specificity"),
    ("calculate_crisproff", "features"),
])
def test_scoring_exports(name: str, module: str) -> None:
    function = getattr(calculator, name)
    assert callable(function)
    assert name in calculator.__all__
    assert function is getattr(import_module(f"DeepOneCalculator.{module}"), name)
    assert function.__name__ == name


@pytest.mark.parametrize(("name", "modules"), [
    ("score_grnas", ("api",)),
    ("calculate_free_energy", ("features", "crispron", "specificity")),
])
def test_old_scoring_names_removed(name: str, modules: tuple[str, ...]) -> None:
    assert name not in calculator.__all__
    assert not hasattr(calculator, name)
    for module in modules:
        assert not hasattr(import_module(f"DeepOneCalculator.{module}"), name)


def test_crisproff_defaults_unchanged() -> None:
    parameters = inspect.signature(calculator.calculate_crisproff).parameters
    assert parameters["off_target_seq"].default is None
    assert parameters["pam_model"].default == "deepone"
    assert parameters["pam_model"].kind is inspect.Parameter.KEYWORD_ONLY
