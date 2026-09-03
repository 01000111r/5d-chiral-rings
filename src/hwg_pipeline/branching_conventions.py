"""Canonical microscopic ``(B,I)`` conventions for branching comparisons.

``B`` is microscopic baryon number, normalized by ``B(Q)=1``.  This module
deliberately never constructs the legacy current-neutral shifted basis.
"""
from dataclasses import dataclass
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path

import yaml
from sage.all import QQ, matrix, vector


class ConventionError(ValueError):
    """A canonical branching convention is invalid or contradicts evidence."""


def exact_rational(value, field="value"):
    if isinstance(value, float):
        raise ConventionError(f"{field} must use exact QQ data, not a float")
    if isinstance(value, dict):
        return QQ(value["numerator"], value["denominator"])
    return QQ(value)


def beta_to_baryon(beta, gauge_rank=3):
    """Convert beta grading to microscopic baryon number for ``SU(Nc)``.

    The default preserves the historical SU(3) API; canonical callers must
    pass the gauge rank supplied by their signed branching specification.
    """
    return exact_rational(gauge_rank, "gauge_rank") * exact_rational(beta, "beta charge")


def negative_cs_instanton_baryon(gauge_rank, signed_k):
    return -exact_rational(gauge_rank) - exact_rational(signed_k, "signed_k")


def conjugate_labels(labels):
    return list(reversed(labels))


def conjugate_operator(operator):
    """Conjugate an operator in a fixed theory; in particular I changes sign."""
    out = dict(operator)
    out["representation"] = conjugate_labels(operator["representation"])
    out["B"], out["I"] = -exact_rational(operator["B"]), -exact_rational(operator["I"])
    return out


def reverse_cs_orientation(operator, new_signed_k, gauge_rank=3):
    """Choose the conjugate *I=+1* current after reversing CS orientation."""
    out = dict(operator)
    out["representation"] = conjugate_labels(operator["representation"])
    out["B"] = negative_cs_instanton_baryon(gauge_rank, new_signed_k)
    out["I"] = QQ.one()
    return out


def render_signed_title(gauge_rank, signed_k, flavours):
    k = exact_rational(signed_k, "signed_k")
    level = "0" if not k else str(k)
    return rf"SU({gauge_rank})_{{{level}}}+{flavours}F"


def render_canonical_charge(representation, B, I):
    """Canonical report fragment; B has no mic/hat/shifted decoration."""
    return rf"[{','.join(map(str, representation))}]_{{B={B},I={I}}}"


def classify_stored_map(instanton_B, instanton_I):
    return "legacy_current_neutral" if exact_rational(instanton_B) == 0 and exact_rational(instanton_I) == 1 else "canonical_candidate"


def solve_two_anchor_map(raw_order, classical, instanton):
    """Solve ``(B,I)=M raw`` directly from two physical anchors over QQ."""
    anchors = (classical, instanton)
    raw = matrix(QQ, [[exact_rational(a["raw_charges"][n]) for a in anchors]
                      for n in raw_order])
    target = matrix(QQ, [[exact_rational(a["target"][n]) for a in anchors]
                         for n in ("B", "I")])
    if raw.rank() != 2:
        raise ConventionError("noninvertible anchor system: anchors do not have rank two")
    solved = target * raw.inverse()
    if any(solved * raw.column(i) != target.column(i) for i in range(2)):
        raise ConventionError("anchor residual does not vanish")
    return solved, solved.inverse()


def _labels(factors):
    return [f["labels"] for f in factors]


def _parent_labels(parent):
    if "parent_factors" in parent:
        return _labels(parent["parent_factors"])
    return parent.get("parent_d5_labels", parent.get("parent_su6_labels"))


def _child_labels(child):
    if "child_factors" in child:
        return child["child_factors"][0]["labels"]
    return child["child_su5_labels"]


def _raw(child):
    if "raw_charges" in child:
        return child["raw_charges"]
    return {"x": child["x_charge"], "q": child["q_charge"]}


def find_anchor(raw_payload, anchor):
    """Representation-aware exact anchor lookup (never raw-charge-only)."""
    matches = []
    for parent in raw_payload["parents"]:
        if parent["degree"] != anchor["degree"] or _parent_labels(parent) != anchor["parent_representations"]:
            continue
        for child in parent["children"]:
            if _child_labels(child) == anchor["child_representation"] and _raw(child) == anchor["raw_charges"]:
                matches.append((parent, child))
    if len(matches) != 1:
        raise ConventionError(f"{anchor['role']}: expected one representation-aware anchor, found {len(matches)}")
    return matches[0]


def _qstr(q):
    q = QQ(q)
    return str(q)


def _matrix_json(m, raw_order):
    return {"raw_charge_order": list(raw_order), "B": [_qstr(x) for x in m.row(0)],
            "I": [_qstr(x) for x in m.row(1)]}


def validate_spec_shape(spec):
    required = {"theory_id", "finite_reference_id", "gauge_rank", "signed_k", "cs_orientation",
                "raw_charge_names", "raw_charge_order", "classical_anchor", "instanton_anchor",
                "expected_charge_map", "charge_lattice", "literature_orientation", "conjugate_check"}
    missing = required - set(spec)
    if missing: raise ConventionError("missing canonical fields: " + ", ".join(sorted(missing)))
    k = exact_rational(spec["signed_k"], "signed_k")
    if k and "absolute_k" in spec: raise ConventionError("unsigned |k| is not canonical; supply signed_k")
    if spec.get("normalization", {}).get("B_of_Q") != 1: raise ConventionError("canonical normalization requires B(Q)=1")
    a = spec["instanton_anchor"]
    if not a.get("child_representation"): raise ConventionError("instanton anchor requires a representation")
    if k and exact_rational(a["target"]["B"]) == 0:
        raise ConventionError("The stored map neutralises an instanton current and is a legacy shifted basis. Supply a signed-k microscopic instanton anchor.")
    if "physical_charge" in spec: raise ConventionError("ambiguous physical_charge field; use explicit B/I keys")
    return k


def _input_paths(root, spec, order):
    uv = root / "generated" / spec["theory_id"] / f"order_{order}"
    fin = root / "generated" / spec["finite_reference_id"] / f"order_{order}"
    raw = uv/"branching_comparison"/"raw_branching.json"
    if not raw.exists(): raw = uv/"manifest_branching"/"branched_refined_plethystic_logarithm.json"
    return [fin/"refined_plethystic_logarithm.json", uv/"refined_plethystic_logarithm.json",
            fin/"reconstruction_checks.json", uv/"reconstruction_checks.json", raw]


def validate_case(root, spec_path, order=10):
    spec = yaml.safe_load(Path(spec_path).read_text())
    k = validate_spec_shape(spec)
    paths = _input_paths(Path(root), spec, order)
    for p in paths:
        if not p.exists(): raise ConventionError(f"missing stored input: {p}")
    for p in (paths[2], paths[3]):
        if not json.loads(p.read_text())["validation_results"]["all_passed"]:
            raise ConventionError(f"failed reconstruction evidence: {p}")
    raw = json.loads(paths[4].read_text())
    classical, instanton = spec["classical_anchor"], spec["instanton_anchor"]
    compact_raw = "coefficients_by_t_degree" in raw
    def compact_find(anchor):
        entries=raw["coefficients_by_t_degree"].get(str(anchor["degree"]),[])
        hits=[e for e in entries if e["child_dynkin_labels"]==anchor["child_representation"]
              and {n:int(e["raw_charges"][n]) for n in spec["raw_charge_order"]}==anchor["raw_charges"]]
        if len(hits)!=1: raise ConventionError(f"{anchor['role']}: expected one exact raw child, found {len(hits)}")
        native=json.loads(paths[1].read_text())["coefficients_by_t_degree"].get(str(anchor["degree"]),[])
        if not any([factor["dynkin_labels"] for factor in e["irreducible_representations"]]==anchor["parent_representations"]
                   if anchor["parent_representations"] and isinstance(anchor["parent_representations"][0], list)
                   else e["irreducible_representations"][0]["dynkin_labels"]==anchor["parent_representations"] for e in native):
            raise ConventionError(f"{anchor['role']}: parent representation absent")
        return hits[0]
    if compact_raw: compact_find(classical); compact_find(instanton)
    else: find_anchor(raw, classical); find_anchor(raw, instanton)
    finite = json.loads(paths[0].read_text())
    beta = exact_rational(classical["finite_beta_charge"])
    if beta_to_baryon(beta, spec["gauge_rank"]) != exact_rational(classical["target"]["B"]):
        raise ConventionError("finite-anchor mismatch: B is not Nc B_beta")
    finite_hits = [e for e in finite["coefficients_by_t_degree"][str(classical["degree"])]
                   if e["irreducible_representations"][0]["dynkin_labels"] == classical["child_representation"]
                   and exact_rational(e["abelian_charges"]["beta"]) == beta]
    if not finite_hits: raise ConventionError("finite-anchor mismatch: beta-graded finite PL term absent")
    M, inverse = solve_two_anchor_map(spec["raw_charge_order"], classical, instanton)
    expected = matrix(QQ, [[exact_rational(x) for x in row] for row in spec["expected_charge_map"]])
    if M != expected: raise ConventionError(f"expected-map mismatch: solved {M}, expected {expected}")
    if k:
        wanted = negative_cs_instanton_baryon(spec["gauge_rank"], k)
        if exact_rational(instanton["target"]["B"]) != wanted or exact_rational(instanton["target"]["I"]) != 1:
            raise ConventionError("zero-mode charge mismatch: require B=-3-k and I=1")
    # The conjugate is checked as a fully identified child, then translated.
    conjugate = compact_find(spec["conjugate_check"]) if compact_raw else find_anchor(raw, spec["conjugate_check"])[1]
    cv = M * vector(QQ, [exact_rational(_raw(conjugate)[n]) for n in spec["raw_charge_order"]])
    if list(cv) != [ -exact_rational(instanton["target"]["B"]), QQ(-1) ]:
        raise ConventionError("conjugation failure: conjugate charges are not (-B,-I)")
    bstep = exact_rational(spec["charge_lattice"]["B_step"])
    count = 0
    children = ([e for entries in raw["coefficients_by_t_degree"].values() for e in entries]
                if compact_raw else [c for parent in raw["parents"] for c in parent["children"]])
    for child in children:
            v = M * vector(QQ, [exact_rational(_raw(child)[n]) for n in spec["raw_charge_order"]])
            if v[0]/bstep not in __import__('sage.all', fromlist=['ZZ']).ZZ or v[1] not in __import__('sage.all', fromlist=['ZZ']).ZZ:
                raise ConventionError(f"charge-lattice failure at raw {_raw(child)}: {(v[0],v[1])}")
            count += 1
    return {"theory_id": spec["theory_id"], "signed_k": _qstr(k), "status": "pass",
            "solved_charge_map": _matrix_json(M, spec["raw_charge_order"]),
            "inverse_charge_map": [[_qstr(x) for x in row] for row in inverse.rows()],
            "instanton_representation": instanton["child_representation"],
            "instanton_B": _qstr(instanton["target"]["B"]), "instanton_I": "1",
            "classical_B": _qstr(classical["target"]["B"]), "classical_I": "0",
            "translated_children": count, "charge_lattice": spec["charge_lattice"],
            "input_paths": [str(p.relative_to(root)) for p in paths] + [str(Path(spec_path).relative_to(root))]}


def run_preflight(root, campaign_path, order=10, strict=True):
    root, campaign_path = Path(root), Path(campaign_path)
    campaign = yaml.safe_load(campaign_path.read_text())
    output = root / campaign["output_directory"]
    output.mkdir(parents=True, exist_ok=True)
    results, failures = [], []
    for item in campaign["cases"]:
        try: results.append(validate_case(root, root/item["specification"], order))
        except Exception as exc:
            failures.append({"theory_id": item["theory_id"], "status":"fail", "diagnostic":str(exc)})
            if strict: break
    payload = {"campaign_id": campaign["campaign_id"], "order": order,
               "canonical_B": "microscopic baryon number normalized by B(Q)=1",
               "all_passed": not failures and len(results)==len(campaign["cases"]), "theories": results, "failures": failures}
    _write_json(output/"preflight_results.json", payload)
    checks = {"all_passed": payload["all_passed"], "case_count":len(results), "failures":failures,
              "no_reports_or_branchings_written":True}
    _write_json(output/"preflight_checks.json", checks)
    hashes = {}
    for case in results:
        for rel in case["input_paths"]: hashes[rel] = sha256((root/rel).read_bytes()).hexdigest()
    _write_json(output/"input_hashes.json", hashes)
    manifest = {"campaign_id":campaign["campaign_id"], "order":order,
                "files":["preflight_results.json","preflight_results.md","preflight_checks.json","preflight_manifest.json","input_hashes.json","command_log.txt"]}
    _write_json(output/"preflight_manifest.json", manifest)
    md=["# Signed branching convention preflight", "", "B is microscopic baryon number normalized by B(Q)=1.", ""]
    md += [f"- **PASS** `{x['theory_id']}`: B coefficients {x['solved_charge_map']['B']}; I coefficients {x['solved_charge_map']['I']}." for x in results]
    md += [f"- **FAIL** `{x['theory_id']}`: {x['diagnostic']}" for x in failures]
    (output/"preflight_results.md").write_text("\n".join(md)+"\n")
    (output/"command_log.txt").write_text(f"./scripts/sage-python -m hwg_pipeline validate-branching-conventions {campaign_path.relative_to(root)} --order {order} --strict\nstatus: {'pass' if payload['all_passed'] else 'fail'}\n")
    if strict and failures: raise ConventionError(failures[0]["diagnostic"])
    return payload


def _write_json(path, data):
    path.write_text(json.dumps(data, indent=2, sort_keys=True)+"\n")
