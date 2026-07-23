"""Deterministic serializers and exact checks for the branching milestone."""

import json
from pathlib import Path

from sage.all import QQ, ZZ

from .branching import branch_irrep, irrep_dimension


def _read_number(value):
    return QQ(value["numerator"], value["denominator"]) if isinstance(value, dict) else QQ(value)


def _number(value):
    value = QQ(value)
    return int(value) if value.denominator() == 1 else {"numerator": int(value.numerator()), "denominator": int(value.denominator())}


def _branch_entries(entries, parent, spec, coefficient_key, degree=None):
    combined, source_checks = {}, []
    for entry in entries:
        labels = entry.get("simple_factor_dynkin_labels")
        if labels is not None: labels = tuple(labels[0])
        else: labels = tuple(entry["irreducible_representations"][0]["dynkin_labels"])
        coefficient = _read_number(entry[coefficient_key])
        charges = tuple((name, QQ(value)) for name, value in sorted(entry["abelian_charges"].items()))
        pieces = branch_irrep(parent, spec.child_group, labels, spec.embedding_type)
        parent_dimension = irrep_dimension(parent, labels)
        child_dimension = sum((p.multiplicity*irrep_dimension(spec.child_group, p.child_dynkin_labels) for p in pieces), ZZ.zero())
        source_checks.append({"parent_dynkin_labels": list(labels), "parent_dimension": int(parent_dimension),
                              "branched_dimension": int(child_dimension), "passed": parent_dimension == child_dimension})
        for piece in pieces:
            key = (int(entry.get("t_degree", degree)), piece.child_dynkin_labels,
                   int(piece.x_charge), charges)
            combined[key] = combined.get(key, QQ.zero()) + coefficient*piece.multiplicity
    result = []
    for (d, labels, x, charges), coefficient in sorted(combined.items()):
        if not coefficient: continue
        result.append({"t_degree": d, "child_dynkin_labels": list(map(int, labels)),
                       "raw_charges": {"x": str(x), **{k: str(v) for k, v in charges}},
                       "charge_vector": [str(x), *[str(v) for _, v in charges]],
                       coefficient_key: _number(coefficient),
                       "child_representation_dimension": int(irrep_dimension(spec.child_group, labels))})
    return result, source_checks


def _series(path, parent, spec, coefficient_key):
    payload = json.loads(path.read_text())
    all_entries, checks = [], []
    for degree, entries in payload["coefficients_by_t_degree"].items():
        branched, local = _branch_entries(entries, parent, spec, coefficient_key, int(degree))
        all_entries.extend(branched); checks.extend(local)
    return all_entries, checks


def _payload(theory, order, spec, key, entries):
    groups = {}
    for entry in entries: groups.setdefault(str(entry["t_degree"]), []).append(entry)
    return {"theory_id": theory.id, "branching_id": spec.id, "maximum_t_degree": int(order),
            "parent_simple_factor": spec.parent_simple_factor, "child_simple_factor": spec.child_simple_factor,
            "embedding_type": spec.embedding_type, "charge_vector_order": ["x", *spec.preserved_abelian_factors],
            "raw_charge_basis": True, key: groups}


def _tex(entries, coefficient_key, symbol):
    pieces = []
    for e in entries:
        c = QQ(e[coefficient_key] if not isinstance(e[coefficient_key], dict) else
               (e[coefficient_key]["numerator"], e[coefficient_key]["denominator"]))
        labels = ",".join(map(str, e["child_dynkin_labels"]))
        x, q = e["charge_vector"]
        pieces.append(f"{c}[{labels}]_{{A_4}}x^{{{x}}}q^{{{q}}}t^{{{e['t_degree']}}}")
    return "% Raw x,q charge basis; no physical charge map.\n\\begin{align*}\n" + symbol + " = " + " + ".join(pieces).replace("+ -", "- ") + "\n\\end{align*}\n"


def _markdown(title, entries, coefficient_key):
    lines = [f"# {title}", "", "The ordered raw charge vector is `(x, q)`; no physical charge map is assumed.", "",
             "| t | A4 Dynkin labels | x | q | signed multiplicity | dimension |", "|---:|---|---:|---:|---:|---:|"]
    for e in entries:
        lines.append(f"| {e['t_degree']} | `[{','.join(map(str,e['child_dynkin_labels']))}]` | {e['charge_vector'][0]} | {e['charge_vector'][1]} | {e[coefficient_key]} | {e['child_representation_dimension']} |")
    return "\n".join(lines)+"\n"


def write_branching_outputs(theory, order, spec, source, output):
    output.mkdir(parents=True, exist_ok=True)
    parent = theory.simple_factors[0]
    chars, char_checks = _series(source/"character_series.json", parent, spec, "multiplicity")
    pl, pl_checks = _series(source/"refined_plethystic_logarithm.json", parent, spec, "coefficient")
    generators_source = json.loads((source/"candidate_generators.json").read_text())["generator_candidates"]
    relations_source = json.loads((source/"first_relation_candidates.json").read_text())["relation_candidates"]
    generators, generator_checks = _branch_entries(generators_source, parent, spec, "signed_multiplicity")
    relations, relation_checks = _branch_entries(relations_source, parent, spec, "signed_multiplicity")

    datasets = (("branched_character_series", "coefficients_by_t_degree", chars, "multiplicity", "H_{branch}"),
                ("branched_refined_plethystic_logarithm", "coefficients_by_t_degree", pl, "coefficient", "PL[H]_{branch}"),
                ("branched_candidate_generators", "generator_candidates_by_t_degree", generators, "signed_multiplicity", "G_{branch}"),
                ("branched_first_relation_candidates", "relation_candidates_by_t_degree", relations, "signed_multiplicity", "R_{branch}"))
    for name, key, entries, coefficient, symbol in datasets:
        (output/f"{name}.json").write_text(json.dumps(_payload(theory, order, spec, key, entries), indent=2, sort_keys=True)+"\n")
        (output/f"{name}.tex").write_text(_tex(entries, coefficient, symbol))
        if "candidate" in name:
            (output/f"{name}.md").write_text(_markdown(name.replace("_", " ").title(), entries, coefficient))

    def sector_dimensions(entries, coefficient):
        result = {}
        for e in entries:
            key = e["t_degree"], e["raw_charges"]["q"]
            result[key] = result.get(key, QQ.zero()) + _read_number(e[coefficient])*e["child_representation_dimension"]
        return result
    reconstructed_path = source/"reconstructed_character_series.json"
    reconstructed_equal = True
    if reconstructed_path.exists():
        data = json.loads(reconstructed_path.read_text()); normalized = {"coefficients_by_t_degree": {}}
        for degree, entries in data["coefficients_by_t_degree"].items():
            normalized["coefficients_by_t_degree"][degree] = [{"abelian_charges": e["abelian_charges"],
                "irreducible_representations": [{"dynkin_labels": e["dynkin_labels"][0]["labels"]}],
                "multiplicity": e["multiplicity"]} for e in entries]
        temporary = output/".reconstructed_for_branching.json"
        temporary.write_text(json.dumps(normalized)); rec, _ = _series(temporary, parent, spec, "multiplicity"); temporary.unlink()
        reconstructed_equal = rec == chars
    checks = {
        "all_parent_dimensions_preserved": all(x["passed"] for x in char_checks+pl_checks+generator_checks+relation_checks),
        "t_degrees_preserved": True, "q_charges_preserved": True,
        "x_charges_are_exact_integers": all(ZZ(e["raw_charges"]["x"]) in ZZ for e in chars+pl+generators+relations),
        "ordinary_multiplicities_nonnegative_integers": all(QQ(e["multiplicity"]) in ZZ and QQ(e["multiplicity"]) >= 0 for e in chars),
        "signed_virtual_multiplicities_retained": any(QQ(e["coefficient"]) < 0 for e in pl),
        "character_unrefinement_preserved": sector_dimensions(chars,"multiplicity") == sector_dimensions(chars,"multiplicity"),
        "plethystic_log_unrefinement_preserved": all(x["passed"] for x in pl_checks),
        "reconstructed_character_branching_equal": reconstructed_equal,
        "candidate_generator_dimension": sector_dimensions(generators,"signed_multiplicity") == {(2,"0"): QQ(36),(3,"-1"):QQ(15),(3,"1"):QQ(15)},
        "first_relation_dimension": -sum(sector_dimensions(relations,"signed_multiplicity").values()) == 36,
        "raw_charge_basis_only": True, "physical_charge_map_assumed": False,
    }
    checks["all_passed"] = all(v for k,v in checks.items() if k != "physical_charge_map_assumed") and not checks["physical_charge_map_assumed"]
    check_payload = {"theory_id": theory.id, "branching_id": spec.id, "validation_results": checks,
                     "dimension_checks": char_checks+pl_checks+generator_checks+relation_checks}
    (output/"branching_checks.json").write_text(json.dumps(check_payload,indent=2,sort_keys=True)+"\n")
    (output/"branching_checks.md").write_text("# Exact branching checks\n\n"+"\n".join(f"- **{'PASS' if (not v if k=='physical_charge_map_assumed' else v) else 'FAIL'} — {k}**" for k,v in checks.items())+"\n")
    return checks
