"""Exact :math:`A_n\to A_{n-1}\times U(1)` character branching.

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
        data["raw_branching_u1_name"], tuple(data["preserved_abelian_factors"]))
    _validate_embedding(spec.parent_simple_factor, spec.child_simple_factor, spec.embedding_type)
    expected = ({"child_dynkin_labels": [1] + [0]*(int(spec.child_rank)-1), "x_charge": 1},
                {"child_dynkin_labels": [0]*int(spec.child_rank), "x_charge": -int(spec.parent_rank)})
    if tuple(data["normalization"]["branches"]) != expected:
        raise ValueError("branching normalization must be fundamental -> fundamental_(+1) + singlet_(-n)")
    return spec


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
    return _branch_cached(parent, labels, child, embedding)


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
