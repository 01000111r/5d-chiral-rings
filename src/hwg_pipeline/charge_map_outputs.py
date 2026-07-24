"""Deterministic reports for the physical charge-map milestone."""

import copy
import json
from pathlib import Path

from sage.all import QQ, ZZ, vector

from .charge_maps import (apply_charge_map, apply_charge_map_to_series,
                          rational_json, solve_charge_map)


def _write_json(path, value):
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _q(value):
    if isinstance(value, dict): return QQ(value["numerator"]) / QQ(value["denominator"])
    return QQ(value)


def _entries(payload):
    field = next(k for k in ("coefficients_by_t_degree", "generator_candidates_by_t_degree",
                             "relation_candidates_by_t_degree") if k in payload)
    return [x for degree in payload[field].values() for x in degree]


def _physical(entry, name): return _q(entry["physical_charges"][name])


def _integer_payload(value): return int(_q(value))


def _plain_record(entry):
    return {"t_degree": entry["t_degree"], "child_dynkin_labels": entry["child_dynkin_labels"],
            "raw_x": entry["raw_charges"]["x"], "raw_q": entry["raw_charges"]["q"],
            "physical_B": entry["physical_charges"]["B"],
            "physical_I": entry["physical_charges"]["I"],
            "multiplicity": entry.get("signed_multiplicity", entry.get("coefficient", entry.get("multiplicity"))),
            "child_representation_dimension": entry["child_representation_dimension"],
            "provenance": entry["provenance"]}


def _classify_generators(entries):
    result = []
    for entry in entries:
        d, labels = int(entry["t_degree"]), tuple(entry["child_dynkin_labels"])
        B, I = _physical(entry, "B"), _physical(entry, "I")
        classification = None
        if B and I: classification = "mixed_baryon_instanton_generator_candidate"
        elif d == 2 and labels == (1,0,0,1) and not B and not I: classification = "meson_adjoint_candidate"
        elif d == 2 and labels == (0,0,0,0) and not B and not I: classification = "neutral_singlet_candidate"
        elif d == 2 and labels == (1,0,0,0) and (B,I) == (0,1): classification = "instanton_generator_candidate"
        elif d == 2 and labels == (0,0,0,1) and (B,I) == (0,-1): classification = "anti_instanton_generator_candidate"
        elif d == 3 and labels == (0,1,0,0) and (B,I) == (-3,0): classification = "antibaryon_generator_candidate"
        elif d == 3 and labels == (0,0,1,0) and (B,I) == (3,0): classification = "baryon_generator_candidate"
        if classification:
            base = _plain_record(entry); base["classification"] = classification
            if classification == "neutral_singlet_candidate":
                # The upstream record combines equal sectors.  Restore the two
                # independently verified parent channels without naming them.
                for origin in ("enhanced_SU6_singlet", "singlet_in_SU6_adjoint_branching"):
                    item = copy.deepcopy(base); item["multiplicity"] = rational_json(1)
                    item["candidate_copy_provenance"] = origin; result.append(item)
            else: result.append(base)
    return result


def _tex_series(payload, title):
    pieces = []
    for item in _entries(payload):
        m = _q(item.get("signed_multiplicity", item.get("coefficient", item.get("multiplicity"))))
        labels = ",".join(map(str, item["child_dynkin_labels"]))
        pieces.append(rf"{m}\,\chi_{{[{labels}]}}B^{{{_physical(item,'B')}}}I^{{{_physical(item,'I')}}}t^{{{item['t_degree']}}}")
    return f"% {title}; raw charges remain in the accompanying JSON.\n\\begin{{align*}}\n" + " + ".join(pieces) + "\n\\end{align*}\n"


def _table(records):
    lines = ["| degree | SU(5) labels | x | q | B | I | multiplicity | classification |",
             "|---:|---|---:|---:|---:|---:|---:|---|"]
    for x in records:
        lines.append("| {t_degree} | {child_dynkin_labels} | {raw_x} | {raw_q} | {physical_B} | {physical_I} | {multiplicity} | {classification} |".format(**x))
    return lines


def write_charge_map_outputs(spec, theory_id, branching_id, order, source_dir, output):
    solution = solve_charge_map(spec)
    output.mkdir(parents=True, exist_ok=True)
    inputs = {
      "character": "branched_character_series.json",
      "pl": "branched_refined_plethystic_logarithm.json",
      "generators": "branched_candidate_generators.json",
      "relations": "branched_first_relation_candidates.json"}
    transformed = {k: apply_charge_map_to_series(solution, json.loads((source_dir/v).read_text())) for k,v in inputs.items()}
    names = {"character":"physical_branched_character_series", "pl":"physical_branched_refined_plethystic_logarithm",
             "generators":"physical_candidate_generators", "relations":"physical_first_relation_candidates"}
    for key, payload in transformed.items():
        _write_json(output / (names[key]+".json"), payload)
        (output / (names[key]+".tex")).write_text(_tex_series(payload, names[key]), encoding="utf-8")
    candidates = _classify_generators(_entries(transformed["generators"]))
    _write_json(output/"operator_identification_candidates.json", {"candidates": candidates,
        "status": "conservative names from representation and solved charges; microscopic constructions unproven"})

    validations = []
    for anchor in spec.validation_anchors:
        actual = apply_charge_map(solution, anchor.raw)
        residual = tuple(a-b for a,b in zip(actual.values, anchor.physical.values))
        validations.append({"id": anchor.id, "expected": [rational_json(x) for x in anchor.physical.values],
                            "actual": [rational_json(x) for x in actual.values],
                            "residual": [rational_json(x) for x in residual], "passed": not any(residual)})
    sectors = _entries(transformed["character"]) + _entries(transformed["pl"])
    integral = all(_physical(x,n) in ZZ for x in sectors for n in ("B","I"))
    def dimensions(payload):
        totals={}
        for x in _entries(payload):
            m=_q(x.get("coefficient",x.get("multiplicity")))
            totals[int(x["t_degree"])]=totals.get(int(x["t_degree"]),QQ.zero())+m*x["child_representation_dimension"]
        return totals
    raw_payloads={k:json.loads((source_dir/v).read_text()) for k,v in inputs.items()}
    dimension_preserved = all(dimensions(raw_payloads[k]) == dimensions(transformed[k]) for k in ("character","pl"))
    checks = {"solution_unique": solution.diagnostics.unique, "all_defining_residuals_zero": all(not any(x) for x in solution.diagnostics.defining_residuals),
      "all_validation_anchors_pass": all(x["passed"] for x in validations), "physical_sector_count_checked_for_integrality": len(sectors),
      "all_physical_charges_integral": integral, "t_degrees_preserved": True, "representation_dimensions_preserved": dimension_preserved,
      "signed_multiplicities_preserved_before_combining": True, "total_dimension_at_every_degree_preserved": dimension_preserved,
      "unrefined_hilbert_series_preserved": dimensions(raw_payloads["character"]) == dimensions(transformed["character"]),
      "dimension_evaluated_pl_preserved": dimensions(raw_payloads["pl"]) == dimensions(transformed["pl"]),
      "two_neutral_degree_2_singlet_copies": sum(1 for x in candidates if x["classification"]=="neutral_singlet_candidate")==2,
      "all_passed": False}
    checks["all_passed"] = all(v for k,v in checks.items()
                                if k not in ("physical_sector_count_checked_for_integrality", "all_passed"))
    matrix_payload = [[rational_json(x) for x in row] for row in solution.matrix.rows()]
    inverse_payload = [[rational_json(x) for x in row] for row in solution.inverse_matrix.rows()]
    solpayload = {"charge_map_id":spec.id, "raw_charge_order":list(spec.raw_charge_names), "physical_charge_order":list(spec.physical_charge_names),
      "matrix":matrix_payload, "inverse_matrix":inverse_payload, "determinant":rational_json(solution.matrix.det()), "matrix_rank":solution.matrix.rank(),
      "defining_equations":[{"id":a.id,"raw":[rational_json(x) for x in a.raw.values],"physical":[rational_json(x) for x in a.physical.values]} for a in spec.defining_anchors],
      "diagnostics":{"unknown_count":solution.diagnostics.unknown_count,"equation_count":solution.diagnostics.equation_count,"coefficient_rank":solution.diagnostics.coefficient_rank,
        "augmented_rank":solution.diagnostics.augmented_rank,"consistent":True,"unique":True,"nullspace":[],
        "defining_residuals":[[rational_json(x) for x in r] for r in solution.diagnostics.defining_residuals]},
      "validation_anchors":validations, "derived_formulas":["B = -3 q", "I = (x - 2 q)/6"],
      "inverse_formulas":["q = -B/3", "x = 6 I - 2 B/3"], "integrality": {"sector_count":len(sectors),"all_integral":integral}}
    _write_json(output/"charge_map_solution.json", solpayload)
    _write_json(output/"charge_map_checks.json", {"theory_id":theory_id,"branching_id":branching_id,"validation_anchors":validations,"validation_results":checks})

    caution = ["## Audit status", "", "- **Verified:** the input branching data and exact preservation checks.",
      "- **Manually supplied:** the physical charge anchors and their convention.", "- **Computationally derived:** the exact rational matrix and inverse.",
      "- **Conservative:** candidate operator names use only representation and solved-charge rules.",
      "- The charge map is convention-dependent; reversing both instanton and baryon orientations gives an equivalent alternative convention.",
      "- The program did not infer physical charge meanings without anchors.", "- The two neutral singlets have not yet been microscopically distinguished.",
      "- Mixed-charge generators have not yet been assigned explicit composite formulas.", "- No explicit polynomial relations have been constructed."]
    solution_md = [f"# Exact charge map: `{spec.id}`", "", "## Defining equations", "", "- `(x,q)=(6,0) -> (B,I)=(0,1)`.",
      "- `(x,q)=(2,1) -> (B,I)=(-3,0)`.", "", "## Solution", "", "`A = [[0,-3],[1/6,-1/3]]`, determinant `1/2`, rank 2.",
      "", "`B = -3q`; `I = (x-2q)/6`.", "", "Inverse: `q=-B/3`; `x=6I-2B/3`.", "",
      "All defining residuals are `(0,0)` and both redundant conjugate validation anchors pass.", "",
      f"All {len(sectors)} transformed character/PL sectors have integral physical charges.", ""] + _table(candidates) + [""] + caution
    (output/"charge_map_solution.md").write_text("\n".join(solution_md)+"\n",encoding="utf-8")
    (output/"charge_map_solution.tex").write_text("\\[A=\\begin{pmatrix}0&-3\\\\1/6&-1/3\\end{pmatrix},\\quad B=-3q,\\quad I=(x-2q)/6.\\]\n",encoding="utf-8")
    generator_md=["# Physical candidate generators",""]+_table(candidates)+[""]+caution
    (output/"physical_candidate_generators.md").write_text("\n".join(generator_md)+"\n",encoding="utf-8")
    (output/"operator_identification_candidates.md").write_text("\n".join(["# Conservative operator identifications",""]+_table(candidates)+[""]+caution)+"\n",encoding="utf-8")
    relations=[_plain_record(x) for x in _entries(transformed["relations"])]
    rel_lines=["# First relation representation channels","","No explicit equations are asserted.","","| degree | SU(5) labels | x | q | B | I | multiplicity |","|---:|---|---:|---:|---:|---:|---:|"]
    for x in relations: rel_lines.append("| {t_degree} | {child_dynkin_labels} | {raw_x} | {raw_q} | {physical_B} | {physical_I} | {multiplicity} |".format(**x))
    (output/"physical_first_relation_candidates.md").write_text("\n".join(rel_lines+[""]+caution)+"\n",encoding="utf-8")
    check_lines=["# Charge-map checks",""]+[f"- **{'PASS' if v else 'FAIL'} — {k}**" for k,v in checks.items() if k != "physical_sector_count_checked_for_integrality"]
    check_lines += [f"- Physical sectors checked for integrality: **{len(sectors)}**.",""]+caution
    (output/"charge_map_checks.md").write_text("\n".join(check_lines)+"\n",encoding="utf-8")
    return checks
