"""Exact D7 -> A6 raw branching benchmarks through t^10."""
import json
from pathlib import Path
from sage.all import QQ
from hwg_pipeline.branching import D_EMBEDDING, branch_irrep, load_branching_spec
from hwg_pipeline.model import SimpleGroupSpec
from hwg_pipeline.sage_backend import irrep_dimension
ROOT=Path(__file__).parents[1]
OUT=ROOT/'generated/su4_nf7_k5o2_infinite/order_10/manifest_branching'
SPEC=load_branching_spec(ROOT/'theories/branchings/su4_nf7_k5o2_to_manifest.yaml')
D7=SimpleGroupSpec('d7','D',7,'SO(14)',tuple(f'mu_{i}' for i in range(1,8)))
A6=SPEC.child_group

def pieces(labels):
 return {(tuple(map(int,p.child_dynkin_labels)),int(p.x_charge)):int(p.multiplicity) for p in branch_irrep(D7,A6,labels,D_EMBEDDING)}

def test_four_exact_d7_benchmarks_and_dimensions():
 cases={
 (0,1,0,0,0,0,0):{((0,1,0,0,0,0),-4):1,((1,0,0,0,0,1),0):1,((0,)*6,0):1,((0,0,0,0,1,0),4):1},
 (0,0,0,0,0,1,0):{((0,0,0,0,0,1),-5):1,((0,0,0,1,0,0),-1):1,((0,1,0,0,0,0),3):1,((0,)*6,7):1},
 (0,0,0,0,0,0,1):{((0,)*6,-7):1,((0,0,0,0,1,0),-3):1,((0,0,1,0,0,0),1):1,((1,0,0,0,0,0),5):1},
 (2,0,0,0,0,0,0):{((2,0,0,0,0,0),-4):1,((1,0,0,0,0,1),0):1,((0,0,0,0,0,2),4):1}}
 for parent,expected in cases.items():
  assert pieces(parent)==expected
  assert irrep_dimension(D7,parent)==sum(m*irrep_dimension(A6,l) for (l,x),m in expected.items())

def key(e): return (tuple(e['child_dynkin_labels']),int(e['raw_charges']['x']),int(e['raw_charges']['q']))
def test_complete_raw_low_degrees_and_all_checks():
 p=json.loads((OUT/'branched_refined_plethystic_logarithm.json').read_text())['coefficients_by_t_degree']
 assert {key(e):QQ(e['coefficient']) for e in p['2']}=={((0,)*6,0,0):2,((1,0,0,0,0,1),0,0):1,((0,1,0,0,0,0),-4,0):1,((0,0,0,0,1,0),4,0):1}
 expected={( (0,0,0,0,0,1),-5,1):1,((0,0,0,1,0,0),-1,1):1,((0,1,0,0,0,0),3,1):1,((0,)*6,7,1):1,((0,)*6,-7,-1):1,((0,0,0,0,1,0),-3,-1):1,((0,0,1,0,0,0),1,-1):1,((1,0,0,0,0,0),5,-1):1,((0,)*6,0,0):-1,((2,0,0,0,0,0),-4,0):-1,((1,0,0,0,0,1),0,0):-1,((0,0,0,0,0,2),4,0):-1}
 assert {key(e):QQ(e['coefficient']) for e in p['4']}==expected
 checks=json.loads((OUT/'branching_checks.json').read_text())['validation_results']; assert checks['all_passed']
 index={(int(d),*key(e)):QQ(e['coefficient']) for d,es in p.items() for e in es}
 assert all(index[(d,tuple(reversed(l)),-x,-q)]==m for (d,l,x,q),m in index.items())
 assert all(e['raw_provenance'] for es in p.values() for e in es)
