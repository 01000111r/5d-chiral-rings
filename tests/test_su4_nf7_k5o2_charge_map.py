"""Canonical physical charges for the complete D7-branched PL."""
import json
from pathlib import Path
import yaml
from sage.all import QQ, ZZ, matrix, vector
from hwg_pipeline.branching_conventions import negative_cs_instanton_baryon, solve_two_anchor_map
from hwg_pipeline.charge_maps import load_charge_map_spec, solve_charge_map
ROOT=Path(__file__).parents[1]
RAW=ROOT/'generated/su4_nf7_k5o2_infinite/order_10/manifest_branching'
PHYS=RAW/'physical_charges'
def q(v): return QQ(v['numerator'])/QQ(v['denominator']) if isinstance(v,dict) else QQ(v)
def entries(name):
 p=json.loads((PHYS/name).read_text()); return [e for es in p['coefficients_by_t_degree'].values() for e in es]
def key(e): return (tuple(e['child_dynkin_labels']),q(e['physical_charges']['B']),q(e['physical_charges']['I']))

def test_two_defining_anchors_exact_solution_inverse_and_zero_modes():
 spec=load_charge_map_spec(ROOT/'theories/charge_maps/su4_nf7_k5o2_manifest_canonical.yaml'); sol=solve_charge_map(spec)
 signed=yaml.safe_load((ROOT/'theories/branching/su4_nf7_k5o2_to_finite.yaml').read_text())
 shared,inv=solve_two_anchor_map(signed['raw_charge_order'],signed['classical_anchor'],signed['instanton_anchor'])
 assert sol.matrix==shared==matrix(QQ,[[QQ(3)/8,QQ(35)/8],[-QQ(1)/4,-QQ(1)/4]])
 assert sol.inverse_matrix==inv==matrix(QQ,[[-QQ(1)/4,-QQ(35)/8],[QQ(1)/4,QQ(3)/8]])
 assert sol.matrix.rank()==2 and sol.matrix.det()==1 and sol.matrix*sol.inverse_matrix==1
 assert negative_cs_instanton_baryon(4,QQ('-5/2'))==QQ('-3/2')
 assert QQ('-3/2')+QQ(7)/2==2 # r=2, hence Lambda^2 7 = [0,1,0,0,0,0]

def test_redundant_anchors_lattice_roundtrip_and_conjugation():
 spec=load_charge_map_spec(ROOT/'theories/charge_maps/su4_nf7_k5o2_manifest_canonical.yaml'); sol=solve_charge_map(spec)
 for a in spec.validation_anchors: assert tuple(sol.matrix*vector(QQ,a.raw.values))==a.physical.values
 es=entries('physical_branched_refined_plethystic_logarithm.json')
 idx={(e['t_degree'],*key(e)):q(e['coefficient']) for e in es}
 for e in es:
  B,I=q(e['physical_charges']['B']),q(e['physical_charges']['I'])
  assert 2*B in ZZ and I in ZZ
  raw=vector(QQ,[q(e['raw_charges'][n]) for n in ('x','q')])
  assert sol.inverse_matrix*vector(QQ,[B,I])==raw
  l=tuple(e['child_dynkin_labels']); assert idx[(e['t_degree'],tuple(reversed(l)),-B,-I)]==q(e['coefficient'])
 assert json.loads((PHYS/'charge_map_checks.json').read_text())['validation_results']['all_passed']

def test_exact_physical_t2_and_t4():
 es=entries('physical_branched_refined_plethystic_logarithm.json'); by={d:{key(e):q(e['coefficient']) for e in es if e['t_degree']==d} for d in (2,4)}
 assert by[2]=={((0,)*6,0,0):2,((1,0,0,0,0,1),0,0):1,((0,1,0,0,0,0),QQ('-3/2'),1):1,((0,0,0,0,1,0),QQ('3/2'),-1):1}
 spec=load_charge_map_spec(ROOT/'theories/charge_maps/su4_nf7_k5o2_manifest_canonical.yaml')
 expected={(tuple(e['dynkin_labels']),q(e['B']),q(e['I'])):q(e['coefficient']) for e in spec.physical_pl_benchmarks[4]}
 assert by[4]==expected
