"""Exact, anchor-defined linear maps between abelian charge bases."""

from dataclasses import dataclass
from pathlib import Path
from typing import Mapping
import copy
import json

import yaml
from sage.all import QQ, matrix, vector


def rational_json(value):
    """Return an exact, unambiguous JSON representation of a rational."""
    value = QQ(value)
    return {"numerator": int(value.numerator()), "denominator": int(value.denominator())}


def _exact(value):
    if isinstance(value, float):
        raise ValueError("floating-point charge data is forbidden")
    if isinstance(value, Mapping):
        return QQ(value["numerator"]) / QQ(value["denominator"])
    return QQ(value)


@dataclass(frozen=True)
class ChargeVector:
    names: tuple[str, ...]
    values: tuple

    def __post_init__(self):
        object.__setattr__(self, "names", tuple(self.names))
        object.__setattr__(self, "values", tuple(_exact(x) for x in self.values))
        if len(self.names) != len(self.values) or len(set(self.names)) != len(self.names):
            raise ValueError("charge names and values must have equal, non-repeated lengths")


@dataclass(frozen=True)
class ChargeAnchor:
    id: str
    raw: ChargeVector
    physical: ChargeVector
    t_degree: int | None = None
    dynkin_labels: tuple[int, ...] = ()
    justification: str = ""


@dataclass(frozen=True)
class ChargeMapSpec:
    id: str
    raw_charge_names: tuple[str, ...]
    physical_charge_names: tuple[str, ...]
    defining_anchors: tuple[ChargeAnchor, ...]
    validation_anchors: tuple[ChargeAnchor, ...] = ()
    notes: tuple[str, ...] = ()
    expected_unique: bool | None = None
    charge_lattice: Mapping | None = None


@dataclass(frozen=True)
class ChargeMapDiagnostics:
    unknown_count: int
    equation_count: int
    coefficient_rank: int
    augmented_rank: int
    consistent: bool
    unique: bool
    nullspace: tuple[tuple, ...]
    defining_residuals: tuple[tuple, ...]


@dataclass(frozen=True)
class ChargeMapSolution:
    spec: ChargeMapSpec
    matrix: object | None
    inverse_matrix: object | None
    diagnostics: ChargeMapDiagnostics


class InconsistentChargeMapError(ValueError):
    def __init__(self, diagnostics):
        super().__init__("charge-map anchors are inconsistent")
        self.diagnostics = diagnostics


def solve_charge_map(spec):
    """Solve ``physical = A raw`` over QQ using defining anchors only."""
    nr, np = len(spec.raw_charge_names), len(spec.physical_charge_names)
    rows, rhs = [], []
    for anchor in spec.defining_anchors:
        if anchor.raw.names != spec.raw_charge_names or anchor.physical.names != spec.physical_charge_names:
            raise ValueError("anchor charge bases do not match specification")
        for i, target in enumerate(anchor.physical.values):
            row = [QQ.zero()] * (nr * np)
            row[i * nr:(i + 1) * nr] = anchor.raw.values
            rows.append(row); rhs.append(target)
    coefficients = matrix(QQ, rows, ncols=nr*np)
    targets = vector(QQ, rhs)
    augmented = coefficients.augment(matrix(QQ, len(rhs), 1, rhs))
    rank, arank = coefficients.rank(), augmented.rank()
    consistent = rank == arank
    unique = consistent and rank == nr*np
    kernel = tuple(tuple(v) for v in coefficients.right_kernel().basis())
    diagnostics = ChargeMapDiagnostics(nr*np, len(rhs), rank, arank, consistent,
                                       unique, kernel, ())
    if not consistent:
        raise InconsistentChargeMapError(diagnostics)
    if not unique:
        return ChargeMapSolution(spec, None, None, diagnostics)
    flat = coefficients.solve_right(targets)
    solved = matrix(QQ, np, nr, flat)
    residuals = tuple(tuple(solved * vector(QQ, a.raw.values) - vector(QQ, a.physical.values))
                      for a in spec.defining_anchors)
    inverse = solved.inverse() if nr == np and solved.is_invertible() else None
    diagnostics = ChargeMapDiagnostics(nr*np, len(rhs), rank, arank, True, True,
                                       kernel, residuals)
    return ChargeMapSolution(spec, solved, inverse, diagnostics)


def apply_charge_map(solution, raw_vector):
    if solution.matrix is None:
        raise ValueError("a unique charge map is required")
    values = raw_vector.values if isinstance(raw_vector, ChargeVector) else tuple(raw_vector)
    if len(values) != len(solution.spec.raw_charge_names):
        raise ValueError("raw charge vector has the wrong dimension")
    return ChargeVector(solution.spec.physical_charge_names,
                        tuple(solution.matrix * vector(QQ, map(_exact, values))))


def apply_charge_map_to_records(solution, records):
    """Transform and combine JSON-like records, retaining raw contributors."""
    groups = {}
    for source in records:
        item = copy.deepcopy(source)
        raw = tuple(_exact(item["raw_charges"][name]) for name in solution.spec.raw_charge_names)
        physical = apply_charge_map(solution, raw)
        item["raw_charges"] = {name: rational_json(value)
                               for name, value in zip(solution.spec.raw_charge_names, raw)}
        item["raw_charge_vector"] = [rational_json(x) for x in raw]
        item["physical_charges"] = {n: rational_json(v) for n, v in zip(physical.names, physical.values)}
        item["physical_charge_vector"] = [rational_json(x) for x in physical.values]
        item.pop("charge_vector", None)
        multkey = "signed_multiplicity" if "signed_multiplicity" in item else ("coefficient" if "coefficient" in item else "multiplicity")
        multiplicity = _exact(item[multkey])
        provenance = copy.deepcopy(item)
        key = (int(item["t_degree"]), tuple(item["child_dynkin_labels"]), physical.values)
        if key not in groups:
            item[multkey] = rational_json(multiplicity)
            item["provenance"] = [provenance]
            groups[key] = [item, multiplicity, multkey]
        else:
            groups[key][1] += multiplicity
            groups[key][0]["provenance"].append(provenance)
    result = []
    for key in sorted(groups):
        item, multiplicity, multkey = groups[key]
        if multiplicity:
            item[multkey] = rational_json(multiplicity)
            item["provenance"].sort(key=lambda x: json.dumps(x, sort_keys=True))
            result.append(item)
    return result


def apply_charge_map_to_series(solution, series):
    if isinstance(series, Mapping):
        output = copy.deepcopy(series)
        field = next(k for k in ("coefficients_by_t_degree", "generator_candidates_by_t_degree",
                                 "relation_candidates_by_t_degree") if k in series)
        transformed = apply_charge_map_to_records(solution,
            [x for degree in series[field].values() for x in degree])
        grouped = {}
        for item in transformed: grouped.setdefault(str(item["t_degree"]), []).append(item)
        output[field] = grouped
        output["raw_charge_basis"] = False
        output["raw_charge_vector_order"] = list(solution.spec.raw_charge_names)
        output["physical_charge_vector_order"] = list(solution.spec.physical_charge_names)
        output.pop("charge_vector_order", None)
        return output
    return apply_charge_map_to_records(solution, series)


def load_charge_map_spec(path):
    data = yaml.safe_load(Path(path).read_text(encoding="utf-8"))
    def anchor(entry):
        return ChargeAnchor(entry["id"], ChargeVector(tuple(data["raw_charges"]), tuple(entry["raw_charges"][x] for x in data["raw_charges"])),
            ChargeVector(tuple(data["physical_charges"]), tuple(entry["physical_charges"][x] for x in data["physical_charges"])),
            entry.get("t_degree"), tuple(entry.get("child_dynkin_labels", entry.get("su5_dynkin_labels", ()))), entry.get("justification", ""))
    return ChargeMapSpec(data["id"], tuple(data["raw_charges"]), tuple(data["physical_charges"]),
        tuple(anchor(x) for x in data["defining_anchors"]), tuple(anchor(x) for x in data["validation_anchors"]),
        tuple(data.get("notes", ())), data.get("expected_unique"), data.get("charge_lattice"))
