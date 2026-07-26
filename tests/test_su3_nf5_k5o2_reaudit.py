"""Understandable exact regression tests for the D5 spinor re-audit."""
import json
from pathlib import Path
from sage.all import QQ, ZZ, WeylCharacterRing, matrix, vector
from hwg_pipeline.branching import branch_irrep, D5_EMBEDDING, validate_d_type_branching
from hwg_pipeline.branching_conventions import find_anchor, solve_two_anchor_map, ConventionError
from hwg_pipeline.model import SimpleGroupSpec
import pytest, yaml
ROOT=Path(__file__).parents[1]
D=SimpleGroupSpec('d','D',5,'SO(10)',tuple('abcde')); A=SimpleGroupSpec('a','A',4,'SU(5)',tuple('abcd'))
def pieces(l): return {(tuple(map(int,z.child_dynkin_labels)),int(z.x_charge)):int(z.multiplicity) for z in branch_irrep(D,A,l,D5_EMBEDDING)}
def test_exact_low_restrictions_and_terminal_distinction():
 assert pieces((0,1,0,0,0))=={((0,1,0,0),-4):1,((1,0,0,1),0):1,((0,0,0,0),0):1,((0,0,1,0),4):1}
 assert pieces((0,0,0,0,1))=={((0,0,0,0),-5):1,((0,0,1,0),-1):1,((1,0,0,0),3):1}
 assert pieces((0,0,0,1,0))=={((0,0,0,1),-3):1,((0,1,0,0),1):1,((0,0,0,0),5):1}
 assert pieces((0,0,0,0,1)) != pieces((0,0,0,1,0))
def test_dimension_would_not_detect_conjugation():
 A4=WeylCharacterRing('A4',style='coroots')
 assert A4((0,1,0,0)).degree()==A4((0,0,1,0)).degree()==10
 assert A4((1,0,0,0)).degree()==A4((0,0,0,1)).degree()==5
 # Nevertheless exact fixed-x labels differ.
 assert ((0,0,1,0),-1) in pieces((0,0,0,0,1)) and ((0,1,0,0),-1) not in pieces((0,0,0,0,1))
 bad=[((0,0,0,0),-5,1),((0,1,0,0),-1,1),((0,0,0,1),3,1)]
 with pytest.raises(ValueError,match='restricted-character failure'):
  validate_d_type_branching('D5',(0,0,0,0,1),bad)
def test_all_reaudited_rules_and_exact_character_points():
 a=json.loads((ROOT/'generated/su3_nf5_k5o2_infinite/order_10/branching_comparison/branching_accuracy_reaudit.json').read_text())
 assert a['rule_count']==78 and a['failure_count']==28
 assert all(r['derived_reconstruction_check'] and all(c['derived_equal'] for c in r['generic_character_checks']) for r in a['rules'])
def test_anchors_map_inverse_conjugation_and_lattice():
 spec=yaml.safe_load((ROOT/'theories/branching/su3_nf5_k5o2_to_finite.yaml').read_text())
 raw=json.loads((ROOT/'generated/su3_nf5_k5o2_infinite/order_10/branching_comparison/raw_branching.json').read_text())
 find_anchor(raw,spec['classical_anchor']); find_anchor(raw,spec['instanton_anchor'])
 old={**spec['classical_anchor'],'parent_representations':[0,0,0,1,0],'raw_charges':{'x':1,'q':1}}
 with pytest.raises(ConventionError): find_anchor(raw,old)
 M,inv=solve_two_anchor_map(spec['raw_charge_order'],spec['classical_anchor'],spec['instanton_anchor'])
 assert M==matrix(QQ,[[QQ(1)/8,-QQ(25)/8],[-QQ(1)/4,QQ(1)/4]])
 assert inv==matrix(QQ,[[-QQ(1)/3,-QQ(25)/6],[-QQ(1)/3,-QQ(1)/6]])
 expected={(-4,0):(-QQ(1)/2,1),(4,0):(QQ(1)/2,-1),(-1,-1):(3,0),(1,1):(-3,0)}
 assert {k:tuple(M*vector(QQ,k)) for k in expected}==expected
 for p in raw['parents']:
  for c in p['children']:
   B,I=M*vector(QQ,[c['x_charge'],c['q_charge']]); assert 2*B in ZZ and I in ZZ
