"""Independent Sage convention and character-restoration tests."""

import json
import pytest
from sage.all import ZZ

from hwg_pipeline import (RepresentationContent, SimpleGroupSpec,
    dimension_refine, expand_hwg, irrep_dimension, load_theory,
    restore_characters, tensor_product, unrefine)


def group(cartan, rank):
    return SimpleGroupSpec("test", cartan, rank, f"{cartan}{rank}",
                           tuple(f"mu_{i}" for i in range(1, rank + 1)))


def test_sage_a2_coroot_convention():
    assert irrep_dimension(group("A", 2), (1, 1)) == ZZ(8)


@pytest.mark.parametrize(("labels", "dimension"), [
    ((0, 0, 0, 0, 0), 1), ((1, 0, 0, 0, 0), 6),
    ((0, 0, 0, 0, 1), 6), ((1, 0, 0, 0, 1), 35),
    ((0, 1, 0, 0, 0), 15), ((0, 0, 0, 1, 0), 15),
    ((0, 1, 0, 1, 0), 189), ((2, 0, 0, 0, 2), 405),
])
def test_sage_a5_dimensions(labels, dimension):
    assert irrep_dimension(group("A", 5), labels) == ZZ(dimension)


def test_a5_fundamental_tensor_antifundamental():
    a5 = group("A", 5)
    decomposition = tensor_product(a5, (1, 0, 0, 0, 0), (0, 0, 0, 0, 1))
    assert decomposition == (((0, 0, 0, 0, 0), ZZ(1)),
                             ((1, 0, 0, 0, 1), ZZ(1)))
    assert 6 * 6 == sum(int(m * irrep_dimension(a5, labels))
                        for labels, m in decomposition) == 1 + 35


def test_representation_content_exact_canonical_arithmetic():
    a5 = group("A", 5); specs = (a5,)
    fundamental = RepresentationContent.single_irrep(specs, ((1, 0, 0, 0, 0),))
    antifund = RepresentationContent.single_irrep(specs, ((0, 0, 0, 0, 1),))
    assert fundamental + fundamental == 2 * fundamental
    assert fundamental * antifund == antifund * fundamental
    assert fundamental - fundamental == RepresentationContent.zero(specs)
    unordered = RepresentationContent(specs, ((((1, 0, 0, 0, 0),), 1),
                                               (((0, 0, 0, 0, 0),), 1)))
    assert list(unordered) == sorted(unordered.terms)


def test_invalid_dynkin_labels_rejected():
    a5 = group("A", 5)
    with pytest.raises(ValueError, match="length"):
        irrep_dimension(a5, (1, 0))
    with pytest.raises(ValueError, match="nonnegative"):
        irrep_dimension(a5, (0, 0, -1, 0, 0))
    with pytest.raises(ValueError, match="floating"):
        irrep_dimension(a5, (0, 0, 1.0, 0, 0))
    with pytest.raises(ValueError, match="floating-point multiplicities"):
        RepresentationContent.single_irrep((a5,), ((0, 0, 0, 0, 0),), 1.0)


def test_physical_character_restoration_and_unrefinement():
    theory = load_theory("theories/su3_nf5_k3o2_infinite.yaml")
    hwg = expand_hwg(theory, 10)
    characters = restore_characters(theory, hwg)
    refined = dimension_refine(characters)
    plain = dict(unrefine(characters))
    assert [plain[ZZ(i)] for i in range(5)] == [1, 0, 36, 30, 630]
    assert all(c in ZZ and c >= 0 for _, c in refined)
    assert plain == {ZZ(d): sum(c for (sd, _), c in refined if sd == d)
                     for d in range(11)}
    assert characters == restore_characters(theory, expand_hwg(theory, 10))


def test_portable_json_shape_has_no_sage_objects():
    content = RepresentationContent.single_irrep((group("A", 2),), ((1, 1),))
    portable = [{"dynkin_labels": [[int(x) for x in labels[0]]],
                 "multiplicity": int(m)} for labels, m in content]
    assert json.loads(json.dumps(portable)) == portable
