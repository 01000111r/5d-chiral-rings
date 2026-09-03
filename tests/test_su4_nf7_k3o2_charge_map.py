import json
from pathlib import Path
from sage.all import QQ, matrix
from hwg_pipeline.charge_maps import load_charge_map_spec, solve_charge_map
from hwg_pipeline.branching_conventions import solve_two_anchor_map
import yaml
ROOT=Path(__file__).parents[1]
def test_su4_charge_map_exact_and_shared():
 s=load_charge_map_spec(ROOT/'theories/charge_maps/su4_nf7_k3o2_manifest_canonical.yaml'); a=solve_charge_map(s)
 c=yaml.safe_load((ROOT/'theories/branching/su4_nf7_k3o2_to_finite.yaml').read_text()); b,_=solve_two_anchor_map(c['raw_charge_order'],c['classical_anchor'],c['instanton_anchor'])
 expected=matrix(QQ,[[-QQ(5)/16,-QQ(49)/16],[QQ(1)/8,-QQ(3)/8]])
 assert a.matrix==b==expected and a.matrix.rank()==2 and a.matrix.det()==QQ(1)/2
 assert a.inverse_matrix==matrix(QQ,[[-QQ(3)/4,QQ(49)/8],[-QQ(1)/4,-QQ(5)/8]])
def test_su4_physical_outputs_pass_all_exact_checks():
 p=json.loads((ROOT/'generated/su4_nf7_k3o2_infinite/order_10/manifest_branching/physical_charges/charge_map_checks.json').read_text())
 assert p['validation_results']['all_passed'] and all(x['passed'] for x in p['validation_anchors'])
