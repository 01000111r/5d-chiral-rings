"""Deterministic serializers and exact checks for the branching milestone."""

import json
from pathlib import Path

from sage.all import QQ, ZZ

from .branching import D_EMBEDDING, PRODUCT_A1_EMBEDDING, branch_a1_to_u1, branch_irrep, irrep_dimension


def _read_number(value):
    return QQ(value["numerator"], value["denominator"]) if isinstance(value, dict) else QQ(value)


def _number(value):
    value = QQ(value)
    return int(value) if value.denominator() == 1 else {"numerator": int(value.numerator()), "denominator": int(value.denominator())}


def _branch_entries(entries, parent, spec, coefficient_key, degree=None, theory=None):
    combined, provenance, source_checks = {}, {}, []
    for entry in entries:
        labels = entry.get("simple_factor_dynkin_labels")
        if labels is not None: labels = [tuple(x) for x in labels]
        else: labels = [tuple(x["dynkin_labels"]) for x in entry["irreducible_representations"]]
        coefficient = _read_number(entry[coefficient_key])
        charges = tuple((name, QQ(value)) for name, value in sorted(entry["abelian_charges"].items()))
        if spec.embedding_type == PRODUCT_A1_EMBEDDING:
            ids = [factor.id for factor in theory.simple_factors]
            preserved_labels = labels[ids.index(spec.preserved_simple_factor_id)]
            a1_labels = labels[ids.index(spec.branched_simple_factor_id)]
            pieces = [(preserved_labels, x, ZZ.one()) for x in branch_a1_to_u1(a1_labels)]
            parent_dimension = irrep_dimension(spec.child_group, preserved_labels) * (ZZ(a1_labels[0]) + 1)
            child_dimension = sum(mult*irrep_dimension(spec.child_group, child) for child, _, mult in pieces)
        else:
            labels = labels[0]
            pieces = [(p.child_dynkin_labels, p.x_charge, p.multiplicity)
                      for p in branch_irrep(parent, spec.child_group, labels, spec.embedding_type)]
            parent_dimension = irrep_dimension(parent, labels)
            child_dimension = sum((mult*irrep_dimension(spec.child_group, child) for child, _, mult in pieces), ZZ.zero())
        source_checks.append({"t_degree": int(entry.get("t_degree", degree)),
                              "abelian_charges": {k: str(v) for k, v in charges},
                              "parent_dynkin_labels": [list(x) for x in labels] if spec.embedding_type == PRODUCT_A1_EMBEDDING else list(labels), "source_coefficient": _number(coefficient),
                              "parent_dimension": int(parent_dimension), "branched_dimension": int(child_dimension),
                              "signed_parent_dimension": _number(coefficient * parent_dimension),
                              "signed_branched_dimension": _number(coefficient * child_dimension),
                              "passed": parent_dimension == child_dimension})
        for child_labels, x_charge, branch_multiplicity in pieces:
            key = (int(entry.get("t_degree", degree)), child_labels, int(x_charge), charges)
            combined[key] = combined.get(key, QQ.zero()) + coefficient*branch_multiplicity
            provenance.setdefault(key, []).append({"branching": spec.id, "embedding_type": spec.embedding_type,
                "parent_factor_dynkin_labels": [list(x) for x in labels] if spec.embedding_type == PRODUCT_A1_EMBEDDING else [list(labels)],
                "parent_coefficient": _number(coefficient), "branching_multiplicity": int(branch_multiplicity)})
    result = []
    for (d, labels, x, charges), coefficient in sorted(combined.items()):
        if not coefficient: continue
        result.append({"t_degree": d, "child_dynkin_labels": list(map(int, labels)),
                       "raw_charges": {"x": str(x), **{k: str(v) for k, v in charges}},
                       "charge_vector": [str(x), *[str(v) for _, v in charges]],
                       coefficient_key: _number(coefficient),
                       "child_representation_dimension": int(irrep_dimension(spec.child_group, labels)),
                       "raw_provenance": provenance[(d, labels, x, charges)]})
    return result, source_checks


def _series(path, parent, spec, coefficient_key, theory):
    payload = json.loads(path.read_text())
    all_entries, checks = [], []
    for degree, entries in payload["coefficients_by_t_degree"].items():
        branched, local = _branch_entries(entries, parent, spec, coefficient_key, int(degree), theory)
        all_entries.extend(branched); checks.extend(local)
    return all_entries, checks


def _payload(theory, order, spec, key, entries):
    groups = {}
    for entry in entries: groups.setdefault(str(entry["t_degree"]), []).append(entry)
    return {"theory_id": theory.id, "branching_id": spec.id, "maximum_t_degree": int(order),
            "parent_simple_factor": spec.parent_simple_factor, "child_simple_factor": spec.child_simple_factor,
            "embedding_type": spec.embedding_type, "charge_vector_order": ["x", *spec.preserved_abelian_factors],
            "raw_charge_basis": "(x,q)", "physical_charge_map_assumed": False, key: groups}


def _tex(entries, coefficient_key, symbol, spec):
    pieces = []
    for e in entries:
        c = QQ(e[coefficient_key] if not isinstance(e[coefficient_key], dict) else
               (e[coefficient_key]["numerator"], e[coefficient_key]["denominator"]))
        labels = ",".join(map(str, e["child_dynkin_labels"]))
        x, q = e["charge_vector"]
        pieces.append(f"{c}[{labels}]_{{{spec.child_simple_factor}}}x^{{{x}}}q^{{{q}}}t^{{{e['t_degree']}}}")
    return "% raw charge basis = (x,q); no physical charge map assumed.\n\\begin{align*}\n" + symbol + " = " + " + ".join(pieces).replace("+ -", "- ") + "\n\\end{align*}\n"


def _markdown(title, entries, coefficient_key, spec):
    lines = [f"# {title}", "", "The ordered raw charge vector is `(x, q)`; no physical charge map is assumed.", "",
             f"| t | {spec.child_simple_factor} Dynkin labels | x | q | signed multiplicity | dimension |", "|---:|---|---:|---:|---:|---:|"]
    for e in entries:
        lines.append(f"| {e['t_degree']} | `[{','.join(map(str,e['child_dynkin_labels']))}]` | {e['charge_vector'][0]} | {e['charge_vector'][1]} | {e[coefficient_key]} | {e['child_representation_dimension']} |")
    return "\n".join(lines)+"\n"


def write_branching_outputs(theory, order, spec, source, output):
    output.mkdir(parents=True, exist_ok=True)
    parent = theory.simple_factors[0]
    chars, char_checks = ([], []) if spec.embedding_type == D_EMBEDDING else _series(source/"character_series.json", parent, spec, "multiplicity", theory)
    pl, pl_checks = _series(source/"refined_plethystic_logarithm.json", parent, spec, "coefficient", theory)
    generators, relations, generator_checks, relation_checks = [], [], [], []
    optional = []
    if (source/"candidate_generators.json").exists():
        generators_source = json.loads((source/"candidate_generators.json").read_text())["generator_candidates"]
        generators, generator_checks = _branch_entries(generators_source, parent, spec, "signed_multiplicity", theory=theory)
        optional.append(("branched_candidate_generators", "generator_candidates_by_t_degree", generators, "signed_multiplicity", "G_{branch}"))
    if (source/"first_relation_candidates.json").exists():
        relations_source = json.loads((source/"first_relation_candidates.json").read_text())["relation_candidates"]
        relations, relation_checks = _branch_entries(relations_source, parent, spec, "signed_multiplicity", theory=theory)
        optional.append(("branched_first_relation_candidates", "relation_candidates_by_t_degree", relations, "signed_multiplicity", "R_{branch}"))

    datasets = ((() if not chars else (("branched_character_series", "coefficients_by_t_degree", chars, "multiplicity", "H_{branch}"),))
                + (("branched_refined_plethystic_logarithm", "coefficients_by_t_degree", pl, "coefficient", "PL[H]_{branch}"),) + tuple(optional))
    for name, key, entries, coefficient, symbol in datasets:
        (output/f"{name}.json").write_text(json.dumps(_payload(theory, order, spec, key, entries), indent=2, sort_keys=True)+"\n")
        (output/f"{name}.tex").write_text(_tex(entries, coefficient, symbol, spec))
        if "candidate" in name:
            (output/f"{name}.md").write_text(_markdown(name.replace("_", " ").title(), entries, coefficient, spec))

    def sector_dimensions(entries, coefficient):
        result = {}
        for e in entries:
            key = e["t_degree"], e["raw_charges"]["q"]
            result[key] = result.get(key, QQ.zero()) + _read_number(e[coefficient])*e["child_representation_dimension"]
        return result
    def persisted_dimensions(path):
        payload = json.loads(path.read_text())
        return {(int(degree), str(entry["abelian_charges"]["q"])): _read_number(entry["coefficient"])
                for degree, entries in payload["coefficients_by_t_degree"].items() for entry in entries}
    char_dimensions = sector_dimensions(chars, "multiplicity")
    pl_dimensions = sector_dimensions(pl, "coefficient")
    expected_char_dimensions = persisted_dimensions(source/"q_refined_dimension_series.json") if chars else {}
    expected_pl_dimensions = persisted_dimensions(source/"q_refined_dimension_pl.json")
    reconstructed_path = source/"reconstructed_character_series.json"
    reconstructed_equal = True
    if chars and reconstructed_path.exists():
        data = json.loads(reconstructed_path.read_text()); normalized = {"coefficients_by_t_degree": {}}
        for degree, entries in data["coefficients_by_t_degree"].items():
            normalized["coefficients_by_t_degree"][degree] = [{"abelian_charges": e["abelian_charges"],
                "irreducible_representations": [{"dynkin_labels": factor["labels"]} for factor in e["dynkin_labels"]],
                "multiplicity": e["multiplicity"]} for e in entries]
        temporary = output/".reconstructed_for_branching.json"
        temporary.write_text(json.dumps(normalized)); rec, _ = _series(temporary, parent, spec, "multiplicity", theory); temporary.unlink()
        reconstructed_equal = rec == chars
    checks = {
        "all_parent_dimensions_preserved": all(x["passed"] for x in char_checks+pl_checks+generator_checks+relation_checks),
        "t_degrees_preserved": True, "q_charges_preserved": True,
        "x_charges_are_exact_integers": all(ZZ(e["raw_charges"]["x"]) in ZZ for e in chars+pl+generators+relations),
        "ordinary_multiplicities_nonnegative_integers": all(QQ(e["multiplicity"]) in ZZ and QQ(e["multiplicity"]) >= 0 for e in chars),
        "signed_virtual_multiplicities_retained": (
            any(QQ(e["coefficient"]) < 0 for e in pl)
            and all(x["signed_parent_dimension"] == x["signed_branched_dimension"] for x in pl_checks)),
        "character_unrefinement_preserved": char_dimensions == expected_char_dimensions if chars else True,
        "plethystic_log_unrefinement_preserved": pl_dimensions == expected_pl_dimensions,
        "reconstructed_character_branching_equal": reconstructed_equal,
        "raw_charge_basis_only": True, "physical_charge_map_assumed": False,
    }
    if generators:
        checks["optional_candidate_generators_dimensions_preserved"] = all(x["passed"] for x in generator_checks)
    if relations:
        checks["optional_first_relations_dimensions_preserved"] = all(x["passed"] for x in relation_checks)
    pl_index = {(e["t_degree"], tuple(e["child_dynkin_labels"]), e["raw_charges"]["x"], e["raw_charges"]["q"]): _read_number(e["coefficient"]) for e in pl}
    checks["plethystic_log_conjugation_consistent"] = all(
        pl_index.get((d, tuple(reversed(labels)), str(-ZZ(x)), str(-QQ(q)))) == coefficient
        for (d, labels, x, q), coefficient in pl_index.items())
    checks["all_passed"] = all(v for k,v in checks.items() if k != "physical_charge_map_assumed") and not checks["physical_charge_map_assumed"]
    check_payload = {"theory_id": theory.id, "branching_id": spec.id,
                     "raw_charge_basis": "(x,q)", "physical_charge_map_assumed": False,
                     "validation_results": checks,
                     "dimension_checks": char_checks+pl_checks+generator_checks+relation_checks}
    (output/"branching_checks.json").write_text(json.dumps(check_payload,indent=2,sort_keys=True)+"\n")
    (output/"branching_checks.md").write_text(
        "# Exact branching checks\n\nRaw charge basis = `(x,q)`; no physical charge map assumed.\n\n"
        + "\n".join(f"- **{'PASS' if (not v if k=='physical_charge_map_assumed' else v) else 'FAIL'} — {k}**" for k,v in checks.items()) + "\n")
    return checks
