"""Cached, validated access to Sage Weyl character rings."""

from functools import lru_cache

from sage.all import CartanType, WeylCharacterRing, ZZ


def _validate_group(group_spec):
    try:
        cartan = CartanType([group_spec.cartan_type, int(group_spec.rank)])
    except (TypeError, ValueError):
        raise ValueError(f"unsupported Cartan type {group_spec.cartan_name}") from None
    if not cartan.is_finite() or not cartan.is_irreducible():
        raise ValueError(f"unsupported Cartan type {group_spec.cartan_name}")
    return group_spec.cartan_name


def _labels(group_spec, dynkin_labels):
    if any(isinstance(x, float) for x in dynkin_labels):
        raise ValueError("floating-point Dynkin labels are forbidden")
    if len(dynkin_labels) != group_spec.rank:
        raise ValueError(f"Dynkin-label length must equal rank {group_spec.rank}")
    labels = []
    for value in dynkin_labels:
        try:
            integer = ZZ(value)
        except (TypeError, ValueError):
            raise ValueError("Dynkin labels must be nonnegative integers") from None
        if integer != value or integer < 0:
            raise ValueError("Dynkin labels must be nonnegative integers")
        labels.append(integer)
    return tuple(labels)


@lru_cache(maxsize=None)
def _ring(cartan_name):
    return WeylCharacterRing(cartan_name, style="coroots")


def character_ring(group_spec):
    return _ring(_validate_group(group_spec))


def irrep(group_spec, dynkin_labels):
    labels = _labels(group_spec, dynkin_labels)
    ring = character_ring(group_spec)
    weight = sum((coefficient * ring.fundamental_weights()[i + 1]
                  for i, coefficient in enumerate(labels)), ring.space().zero())
    return ring(weight)


def irrep_dimension(group_spec, dynkin_labels):
    return ZZ(irrep(group_spec, dynkin_labels).degree())


def _weight_labels(ring, weight):
    coroots = weight.parent().simple_coroots()
    return tuple(ZZ(weight.scalar(coroots[i])) for i in range(1, ring.cartan_type().rank() + 1))


def tensor_product(group_spec, left_labels, right_labels):
    """Return the canonical ``(Dynkin labels, multiplicity)`` decomposition."""
    ring = character_ring(group_spec)
    product = irrep(group_spec, left_labels) * irrep(group_spec, right_labels)
    terms = ((_weight_labels(ring, weight), ZZ(mult))
             for weight, mult in product.monomial_coefficients().items())
    return tuple(sorted(terms))


def symmetric_power(group_spec, dynkin_labels, n):
    """Decompose an irreducible symmetric power using Sage characters."""
    labels = _labels(group_spec, dynkin_labels)
    n = ZZ(n)
    if n <= 0:
        raise ValueError("symmetric-power exponent must be positive")
    character = irrep(group_spec, labels).symmetric_power(int(n))
    ring = character_ring(group_spec)
    return tuple(sorted((_weight_labels(ring, weight), ZZ(mult))
                        for weight, mult in character.monomial_coefficients().items()))


@lru_cache(maxsize=None)
def _adams_decomposition(cartan_name, dynkin_labels, k):
    """Cached irreducible decomposition of a character Adams operation."""
    ring = _ring(cartan_name)
    weight = sum((coefficient * ring.fundamental_weights()[i + 1]
                  for i, coefficient in enumerate(dynkin_labels)), ring.space().zero())
    character = ring(weight).adams_operator(k)
    return tuple(sorted((_weight_labels(ring, weight), ZZ(mult))
                        for weight, mult in character.monomial_coefficients().items()))


def adams_decomposition(group_spec, dynkin_labels, k):
    labels = _labels(group_spec, dynkin_labels)
    k = ZZ(k)
    if k <= 0:
        raise ValueError("Adams index must be a positive integer")
    return _adams_decomposition(_validate_group(group_spec), labels, int(k))
