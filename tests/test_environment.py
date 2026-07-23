"""Checks for the SageMath execution environment and exact arithmetic."""

from sage.all import QQ, WeylCharacterRing


def test_sage_imports():
    """Sage's public interface is available to the test process."""
    import sage.all

    assert sage.all is not None


def test_a2_adjoint_representation_has_dimension_eight():
    """The A2 representation with highest weight [1, 1] is the adjoint."""
    a2 = WeylCharacterRing("A2", style="coroots")

    assert a2(1, 1).degree() == 8


def test_rational_field_parses_three_halves_exactly():
    """A rational string remains an exact element of QQ."""
    value = QQ("3/2")

    assert value == QQ(3) / QQ(2)
    assert value.parent() is QQ
