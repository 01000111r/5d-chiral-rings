"""Deterministic reports for the conservative PL analysis milestone."""

import json

from sage.all import QQ, ZZ

from .operators import (candidate_generators, enumerate_quadratic_channels,
                        extract_operator_content, first_negative_degree,
                        first_relation_candidates)
from .sage_backend import irrep_dimension


def _rat(value):
    return str(value)


def _record(record):
    return {"t_degree": int(record.t_degree),
            "abelian_charges": {k: _rat(v) for k, v in record.abelian_charges},
            "simple_factor_dynkin_labels": [list(map(int, x)) for x in record.dynkin_labels],
            "signed_multiplicity": int(record.signed_multiplicity),
            "representation_dimension": int(record.representation_dimension),
            "sign": record.sign, "classification": record.classification,
            "source_theory_id": record.source_theory_id,
            "mixed_degree": record.mixed_degree}


def _term(labels, multiplicity, theory):
    return {"simple_factor_dynkin_labels": [list(map(int, x)) for x in labels],
            "multiplicity": int(multiplicity),
            "representation_dimension": int(ZZ.prod(irrep_dimension(f, x)
                for f, x in zip(theory.simple_factors, labels)))}


def _coefficient(entries):
    result = {}
    for labels, multiplicity in entries:
        result[labels] = result.get(labels, ZZ(0)) + ZZ(multiplicity)
    return {k: v for k, v in result.items() if v}


def _actual_degree(theory, character_payload, degree, charges):
    result = {}
    for entry in character_payload["coefficients_by_t_degree"].get(str(degree), []):
        q = tuple((f.id, QQ(entry["abelian_charges"][f.id]))
                  for f in theory.abelian_factors)
        if q != charges:
            continue
        by_id = {x["cartan_factor_id"]: tuple(x["dynkin_labels"])
                 for x in entry["irreducible_representations"]}
        labels = tuple(by_id[f.id] for f in theory.simple_factors)
        result[labels] = result.get(labels, ZZ(0)) + ZZ(entry["multiplicity"])
    return result


def write_operator_outputs(theory, order, pl, character_payload, output):
    records = extract_operator_content(pl); generators = candidate_generators(pl)
    relations = first_relation_candidates(pl); degree = first_negative_degree(pl)
    channels = enumerate_quadratic_channels(theory, generators, relations)
    charges = tuple((x.id, ZZ(0)) for x in theory.abelian_factors)
    relation_charges = {x.abelian_charges for x in relations}
    if len(relation_charges) == 1:
        charges = next(iter(relation_charges))
    free = _coefficient((labels, mult) for channel in channels
                        for labels, mult in channel["decomposition"])
    actual = _actual_degree(theory, character_payload, degree, charges)
    deficit = {key: free.get(key, 0) - actual.get(key, 0) for key in set(free) | set(actual)}
    deficit = {k: v for k, v in deficit.items() if v}
    negative = {x.dynkin_labels: -x.signed_multiplicity for x in relations
                if x.abelian_charges == charges}
    channel_json = []
    for channel in channels:
        channel_json.append({"left_generator": _record(channel["left"]),
          "right_generator": _record(channel["right"]), "product_type": channel["product_type"],
          "total_t_degree": int(channel["t_degree"]),
          "total_abelian_charge": {k: _rat(v) for k,v in channel["abelian_charges"]},
          "complete_irreducible_decomposition": [_term(x,m,theory) for x,m in channel["decomposition"]],
          "first_relation_multiplicities": [_term(x,m,theory) for x,m in channel["relation_multiplicities"]]})
    base = {"theory_id": theory.id, "maximum_t_degree": int(order)}
    def write(name, value):
        (output/name).write_text(json.dumps(value, indent=2, sort_keys=True)+"\n")
    write("operator_content.json", {**base, "first_negative_degree": int(degree), "operators": list(map(_record, records))})
    write("candidate_generators.json", {**base, "generator_candidates": list(map(_record, generators))})
    write("first_relation_candidates.json", {**base, "first_negative_degree": int(degree), "relation_candidates": list(map(_record, relations))})
    comparison = {"free_symmetric_algebra_coefficient": [_term(k,v,theory) for k,v in sorted(free.items())],
                  "actual_hilbert_series_coefficient": [_term(k,v,theory) for k,v in sorted(actual.items())],
                  "exact_deficit": [_term(k,v,theory) for k,v in sorted(deficit.items())]}
    write("first_relation_channels.json", {**base, "first_relation_degree": int(degree),
          "quadratic_channels": channel_json, **comparison})
    checks = {"first_negative_degree_exists": degree is not None,
              "deficit_equals_negative_pl_content": deficit == negative,
              "relation_degree_not_mixed": not any(x.mixed_degree for x in relations),
              "dimensions_exact": all(x.representation_dimension in ZZ for x in records)}
    checks["all_passed"] = all(checks.values())
    write("operator_content_checks.json", {**base, "validation_results": checks})
    rep = lambda labels: " ".join("["+",".join(map(str,x))+"]" for x in labels)
    texrow = lambda x: f"{x.signed_multiplicity}\\,{rep(x.dynkin_labels)}\\,t^{{{x.t_degree}}}"
    (output/"operator_content.tex").write_text("% Conservatively classified signed PL content.\n\\begin{align*}\n"+" \\\\\n".join(texrow(x) for x in records)+"\n\\end{align*}\n")
    (output/"candidate_generators.tex").write_text("% Low-degree generator candidates only.\n\\begin{align*}\n"+" \\\\\n".join(texrow(x) for x in generators)+"\n\\end{align*}\n")
    (output/"first_relation_candidates.tex").write_text("% First-relation candidates only.\n\\begin{align*}\n"+" \\\\\n".join(texrow(x) for x in relations)+"\n\\end{align*}\n")
    md = [f"# Operator content: `{theory.id}`", "", "## 1. Verified series data", "",
          f"The first negative refined-PL degree is **{degree}**. Exact reconstruction checks are an input prerequisite.", "",
          "## 2. Computational classification", "",
          "Positive terms below that degree are low-degree generator candidates; negative terms at it are first-relation candidates.", "",
          "## 3. Cautious physical interpretation", "",
          "These labels are computational classifications, not a proof of a minimal coordinate-ring presentation. Later PL terms have not been proven to be minimal generators, defining relations, or syzygies.", ""]
    (output/"operator_content.md").write_text("\n".join(md)+"\n")
    cmd = ["# First relation channels", "", f"Relation degree: **{degree}**.", ""]
    for x in channel_json:
        pieces = " + ".join((str(y["multiplicity"])+"*" if y["multiplicity"] != 1 else "")+rep(tuple(tuple(z) for z in y["simple_factor_dynkin_labels"])) for y in x["complete_irreducible_decomposition"])
        cmd.append(f"- **{x['product_type']}**: {pieces}")
    cmd += ["", "The free coefficient, actual Hilbert coefficient, and their exact deficit are recorded without constructing polynomial equations.", ""]
    (output/"first_relation_channels.md").write_text("\n".join(cmd))
    (output/"operator_content_checks.md").write_text("# Operator-content checks\n\n"+"\n".join(f"- **{'PASS' if v else 'FAIL'} — {k}**" for k,v in checks.items())+"\n")
    return checks
