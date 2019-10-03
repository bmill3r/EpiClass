import pytest


@pytest.mark.parametrize("a", [1])
def test_mutated_base_ambiguous(a):
    assert a > 0
