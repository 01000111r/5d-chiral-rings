"""Exact character branching with a retained raw :math:`U(1)` charge.

The implementation restricts every Sage weight.  In ambient A_n coordinates
the U(1) cocharacter is ``diag(1,...,1,-n)``; no physical interpretation is
attached to this raw charge.
"""

from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

import yaml
from itertools import product

from sage.all import QQ, ZZ

from .model import SimpleGroupSpec
from .sage_backend import irrep, irrep_dimension


EMBEDDING = "A_n_to_A_n_minus_1_u1"
PRODUCT_A1_EMBEDDING = "product_preserve_A_branch_A1_to_u1"
D5_EMBEDDING = "D5_to_A4_u1"
D6_EMBEDDING = "D6_to_A5_u1"


@dataclass(frozen=True)
class BranchingSpec:
    id: str
    source_theory_id: str
    parent_simple_factor: str
    child_simple_factor: str
    embedding_type: str
    child_display_name: str
    raw_branching_u1_name: str
    preserved_abelian_factors: tuple[str, ...]
    preserved_simple_factor_id: str | None = None
    branched_simple_factor_id: str | None = None

    @property
    def parent_rank(self): return ZZ(self.parent_simple_factor[1:])
    @property
    def child_rank(self): return ZZ(self.child_simple_factor[1:])
    @property
    def child_group(self):
        return SimpleGroupSpec("manifest", "A", self.child_rank,
                               self.child_display_name,
                               tuple(f"nu_{i}" for i in range(1, int(self.child_rank)+1)))


@dataclass(frozen=True)
class BranchedIrrepTerm:
    child_dynkin_labels: tuple
    x_charge: object
    multiplicity: object

    def __post_init__(self):
        object.__setattr__(self, "child_dynkin_labels", tuple(ZZ(x) for x in self.child_dynkin_labels))
        object.__setattr__(self, "x_charge", ZZ(self.x_charge))
        object.__setattr__(self, "multiplicity", ZZ(self.multiplicity))


@dataclass(frozen=True)
class BranchedContentTerm:
    t_degree: object
    child_dynkin_labels: tuple
    x_charge: object
    abelian_charges: tuple
    multiplicity: object

    def __post_init__(self):
        object.__setattr__(self, "t_degree", ZZ(self.t_degree))
        object.__setattr__(self, "child_dynkin_labels", tuple(ZZ(x) for x in self.child_dynkin_labels))
        object.__setattr__(self, "x_charge", ZZ(self.x_charge))
        object.__setattr__(self, "abelian_charges", tuple((str(k), QQ(v)) for k, v in self.abelian_charges))
        object.__setattr__(self, "multiplicity", QQ(self.multiplicity))


def load_branching_spec(path):
    data = yaml.safe_load(Path(path).read_text(encoding="utf-8"))
    spec = BranchingSpec(data["id"], data["source_theory_id"], data["parent_simple_factor"],
        data["child_simple_factor"], data["embedding_type"], data["child_display_name"],
        data["raw_branching_u1_name"], tuple(data["preserved_abelian_factors"]),
        data.get("preserved_simple_factor_id"), data.get("branched_simple_factor_id"))
    if spec.embedding_type == PRODUCT_A1_EMBEDDING:
        if spec.parent_simple_factor != "A1" or not spec.child_simple_factor.startswith("A"):
            raise ValueError("product A1 branching requires an A1 branched factor and preserved type-A child")
        if not spec.preserved_simple_factor_id or not spec.branched_simple_factor_id:
            raise ValueError("product branching must identify preserved and branched factor ids")
        if data.get("normalization") != {"su2_irrep_weights": "m,m-2,...,-m"}:
            raise ValueError("A1 -> U(1) normalization must be the exact SU(2) weight convention")
        return spec
    _validate_embedding(spec.parent_simple_factor, spec.child_simple_factor, spec.embedding_type)
    expected = ({"child_dynkin_labels": [1] + [0]*(int(spec.child_rank)-1), "x_charge": 1},
                {"child_dynkin_labels": [0]*int(spec.child_rank), "x_charge": -int(spec.parent_rank)})
    expected_parent = [1] + [0] * (int(spec.parent_rank) - 1)
    if data["normalization"].get("parent_fundamental") != expected_parent:
        raise ValueError("branching normalization must identify the parent fundamental")
    if tuple(data["normalization"]["branches"]) != expected:
        raise ValueError("branching normalization must be fundamental -> fundamental_(+1) + singlet_(-n)")
    return spec


def branch_a1_to_u1(labels):
    """Return the exact weights of the SU(2) irrep with Dynkin label ``m``."""
    if len(labels) != 1 or ZZ(labels[0]) < 0:
        raise ValueError("A1 branching requires one nonnegative Dynkin label")
    m = ZZ(labels[0])
    return tuple(ZZ(x) for x in range(int(m), -int(m)-1, -2))


def _cartan_name(group):
    if isinstance(group, str): return group.replace("_", "")
    return group.cartan_name


def _validate_embedding(parent, child, embedding):
    if embedding != EMBEDDING: raise ValueError(f"unsupported embedding {embedding!r}")
    if not parent.startswith("A") or not child.startswith("A"):
        raise ValueError("embedding requires type A parent and child")
    try: n, m = int(parent[1:]), int(child[1:])
    except ValueError: raise ValueError("invalid Cartan type") from None
    if n < 2 or m != n-1: raise ValueError("embedding requires A_n -> A_(n-1) with n >= 2")
    return n


def _weight_labels(ring, weight):
    roots = weight.parent().simple_coroots()
    return tuple(ZZ(weight.scalar(roots[i])) for i in range(1, ring.cartan_type().rank()+1))


@lru_cache(maxsize=None)
def _branch_cached(parent_name, labels, child_name, embedding):
    n = _validate_embedding(parent_name, child_name, embedding)
    if len(labels) != n or any(ZZ(x) < 0 for x in labels):
        raise ValueError(f"Dynkin-label length must equal parent rank {n}")
    # Validate the representations eagerly through Sage's Weyl machinery.
    parent_spec = SimpleGroupSpec("parent", "A", n, "", tuple(f"p{i}" for i in range(n)))
    child_spec = SimpleGroupSpec("child", "A", n-1, "", tuple(f"c{i}" for i in range(n-1)))
    irrep(parent_spec, labels)

    # A_n highest weights are partitions lambda (last part zero).  Restriction
    # from GL(n+1) to GL(n) is the exact Gelfand--Tsetlin interlacing rule
    # lambda_i >= mu_i >= lambda_(i+1), with multiplicity one.  Passing from
    # mu to its differences gives the A_(n-1) highest weight.  This is the
    # closed-form reconstruction of the same restricted Sage weight system.
    lam = tuple(sum(labels[i:], ZZ.zero()) for i in range(n)) + (ZZ.zero(),)
    boxes = sum(lam, ZZ.zero())
    result = []
    ranges = tuple(range(int(lam[i+1]), int(lam[i])+1) for i in range(n))
    for mu_raw in product(*ranges):
        mu = tuple(ZZ(x) for x in mu_raw)
        child_labels = tuple(mu[i]-mu[i+1] for i in range(n-1))
        charge = (n+1)*sum(mu, ZZ.zero()) - n*boxes
        irrep(child_spec, child_labels)
        result.append(BranchedIrrepTerm(child_labels, charge, 1))
    return tuple(sorted(result, key=lambda term: (term.x_charge, term.child_dynkin_labels)))


def branch_irrep(parent_group, child_group, labels, embedding=EMBEDDING):
    """Return the finite exact irreducible branching of one parent irrep."""
    parent, child = _cartan_name(parent_group), _cartan_name(child_group)
    labels = tuple(ZZ(x) for x in labels)
    if embedding == D5_EMBEDDING:
        if parent != "D5" or child != "A4":
            raise ValueError("D5 embedding requires D5 -> A4")
        return _branch_d5_a4(labels)
    if embedding == D6_EMBEDDING:
        if parent != "D6" or child != "A5":
            raise ValueError("D6 embedding requires D6 -> A5")
        return _branch_d6_a5(labels)
    return _branch_cached(parent, labels, child, embedding)


def validate_d_type_branching(parent_name, labels, actual_children):
    """Strict exact validation; dimensions alone are deliberately insufficient.

    ``actual_children`` is an iterable of ``(A-labels, x, multiplicity)``.
    The expected label/charge multiset is obtained from exact restricted
    weights, not from dimensions.  The resulting error identifies both sides.
    """
    rank=int(str(parent_name)[1:])
    if parent_name not in ("D5", "D6"):
        raise ValueError("strict D-type validation supports the configured D5/D6 embeddings")
    child_rank=rank-1
    parent=SimpleGroupSpec("strict-parent","D",rank,f"SO({2*rank})",tuple(f"p{i}" for i in range(rank)))
    child=SimpleGroupSpec("strict-child","A",child_rank,f"SU({rank})",tuple(f"c{i}" for i in range(child_rank)))
    embedding=D5_EMBEDDING if rank==5 else D6_EMBEDDING
    expected=branch_irrep(parent,child,tuple(labels),embedding)
    e=sorted((tuple(z.child_dynkin_labels),ZZ(z.x_charge),ZZ(z.multiplicity)) for z in expected)
    a=sorted((tuple(ZZ(v) for v in lab),ZZ(x),ZZ(m)) for lab,x,m in actual_children)
    pd=irrep_dimension(parent,labels); ad=sum(m*irrep_dimension(child,lab) for lab,x,m in a)
    if pd != ad or e != a:
        raise ValueError(f"strict restricted-character failure for {parent_name} parent {tuple(labels)}: "
                         f"expected children {e}; actual children {a}; dimension {pd} versus {ad}; "
                         "fixed-charge restricted-weight evidence differs")
    return {"parent":tuple(labels),"expected_children":tuple(e),"actual_children":tuple(a),
            "dimension_equal":True,"restricted_character_equal":True}


@lru_cache(maxsize=None)
def _branch_d5_a4(labels):
    """Restrict every D5 weight and exactly decompose each charge slice.

    In orthonormal D5 coordinates the embedded A4 has the usual coordinate
    differences and ``x=-2*sum(weight)``.  Thus the vector restricts as
    ``5_(-2) + anti-5_(+2)``; this also fixes both spinor-node conventions.
    """
    from sage.all import WeylCharacterRing
    if len(labels) != 5 or any(x < 0 for x in labels):
        raise ValueError("Dynkin-label length must equal parent rank 5")
    D = WeylCharacterRing("D5", style="coroots")
    A = WeylCharacterRing("A4", style="coroots")
    def character(ring, dynkin):
        weight = sum((a*ring.fundamental_weights()[i+1]
                      for i,a in enumerate(dynkin)), ring.space().zero())
        return ring(weight)
    slices = {}
    # ``style="coroots"`` already uses the public Bourbaki D5 ordering.  In
    # particular its fourth and fifth fundamental weights are the two
    # distinct half-spinors.  No terminal-node swap or charge-dependent A4
    # conjugation is permissible at this boundary.
    for weight, multiplicity in character(D, labels).weight_multiplicities().items():
        coords = tuple(weight[i] for i in range(5))
        x = ZZ(-2 * sum(coords))
        dynkin_weight = tuple(ZZ(coords[i] - coords[i+1]) for i in range(4))
        slices.setdefault(x, {})[dynkin_weight] = ZZ(multiplicity)

    answer = []
    for x, remaining in slices.items():
        remaining = {k:v for k,v in remaining.items() if v}
        while remaining:
            dominant = [k for k,v in remaining.items() if v > 0 and all(a >= 0 for a in k)]
            if not dominant:
                raise ValueError("exact A4 character decomposition failed")
            # Pairing with rho is strictly increasing in the dominance order.
            highest = max(dominant, key=lambda k: (sum((i+1)*(5-i)*k[i] for i in range(4)), k))
            coefficient = remaining[highest]
            answer.append(BranchedIrrepTerm(highest, x, coefficient))
            for weight, mult in character(A, highest).weight_multiplicities().items():
                c = tuple(weight[i] for i in range(5))
                key = tuple(ZZ(c[i]-c[i+1]) for i in range(4))
                remaining[key] = remaining.get(key, ZZ.zero()) - coefficient*ZZ(mult)
                if not remaining[key]: del remaining[key]
            if any(v < 0 for v in remaining.values()):
                raise ValueError("non-character encountered in exact A4 decomposition")
    return tuple(sorted(answer, key=lambda term: (term.x_charge, term.child_dynkin_labels)))


@lru_cache(maxsize=None)
def _branch_d6_a5(labels):
    """Exact ``SO(12) -> SU(6) x U(1)`` restriction.

    The project Dynkin-node convention fixes node six by
    ``32' -> 6_(-2) + 20_0 + anti-6_(+2)``.  In orthonormal D6
    coordinates this is ``x=-sum(weight)``.  Charge slices are decomposed
    into complete A5 characters, so this remains exact for every irrep used
    by the order-ten plethystic logarithm.
    """
    from sage.all import WeylCharacterRing
    if len(labels) != 6 or any(x < 0 for x in labels):
        raise ValueError("Dynkin-label length must equal parent rank 6")
    D = WeylCharacterRing("D6", style="coroots")
    A = WeylCharacterRing("A5", style="coroots")
    def character(ring, dynkin):
        weight = sum((a*ring.fundamental_weights()[i+1]
                      for i,a in enumerate(dynkin)), ring.space().zero())
        return ring(weight)
    # Match the repository's public ordering of the two D6 spinor nodes.
    sage_labels = labels[:4] + (labels[5], labels[4])
    slices = {}
    for weight, multiplicity in character(D, sage_labels).weight_multiplicities().items():
        coords = tuple(weight[i] for i in range(6))
        x = ZZ(-sum(coords))
        dynkin_weight = tuple(ZZ(coords[i] - coords[i+1]) for i in range(5))
        slices.setdefault(x, {})[dynkin_weight] = ZZ(multiplicity)
    answer = []
    for x, remaining in slices.items():
        remaining = {k:v for k,v in remaining.items() if v}
        while remaining:
            dominant = [k for k,v in remaining.items() if v > 0 and all(a >= 0 for a in k)]
            if not dominant:
                raise ValueError("exact A5 character decomposition failed")
            highest = max(dominant, key=lambda k: (sum((i+1)*(6-i)*k[i] for i in range(5)), k))
            coefficient = remaining[highest]
            # The adjoint fixes the public A5 orientation opposite to the
            # ambient coordinate-difference labels.  This is one orientation
            # for every weight coset: making spinorial cosets an exception
            # conjugates 6 and anti-6 at fixed x while preserving dimensions.
            child = tuple(reversed(highest))
            answer.append(BranchedIrrepTerm(child, x, coefficient))
            for weight, mult in character(A, highest).weight_multiplicities().items():
                c = tuple(weight[i] for i in range(6))
                key = tuple(ZZ(c[i]-c[i+1]) for i in range(5))
                remaining[key] = remaining.get(key, ZZ.zero()) - coefficient*ZZ(mult)
                if not remaining[key]: del remaining[key]
            if any(v < 0 for v in remaining.values()):
                raise ValueError("non-character encountered in exact A5 decomposition")
    return tuple(sorted(answer, key=lambda term: (term.x_charge, term.child_dynkin_labels)))


def branch_representation_content(content, branching_spec, t_degree=0, abelian_charges=()):
    """Branch linear representation content, preserving every external charge."""
    if len(content.group_specs) != 1: raise ValueError("this milestone branches one parent simple factor")
    terms = []
    for (labels,), coefficient in content:
        for piece in branch_irrep(content.group_specs[0], branching_spec.child_group, labels,
                                  branching_spec.embedding_type):
            terms.append(BranchedContentTerm(t_degree, piece.child_dynkin_labels, piece.x_charge,
                                             abelian_charges, coefficient*piece.multiplicity))
    return _combine_terms(terms)


def _combine_terms(terms):
    combined = {}
    for term in terms:
        key = term.t_degree, term.child_dynkin_labels, term.x_charge, term.abelian_charges
        combined[key] = combined.get(key, QQ.zero()) + term.multiplicity
    return tuple(BranchedContentTerm(*key, coefficient) for key, coefficient in sorted(combined.items()) if coefficient)


def branch_character_series(series, branching_spec):
    terms = []
    for (degree, charges), content in series:
        terms.extend(branch_representation_content(content, branching_spec, degree, charges))
    result = _combine_terms(terms)
    if any(term.multiplicity not in ZZ or term.multiplicity < 0 for term in result):
        raise ValueError("ordinary branched character multiplicities must be nonnegative integers")
    return result


def branch_virtual_series(series, branching_spec):
    terms = []
    for (degree, charges), content in series:
        terms.extend(branch_representation_content(content, branching_spec, degree, charges))
    return _combine_terms(terms)


def branched_dimension(terms, child_group):
    return sum((term.multiplicity*irrep_dimension(child_group, term.child_dynkin_labels)
                for term in terms), QQ.zero())
