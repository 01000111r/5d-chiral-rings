"""Focused mathematical tests for canonical microscopic branching charges."""
from copy import deepcopy
import json
from pathlib import Path
import pytest
import yaml
from sage.all import QQ, matrix
from hwg_pipeline.branching_conventions import *

ROOT=Path(__file__).parents[1]
CAMPAIGN=ROOT/'theories/campaigns/su3_nf5_nf6_signed_branching.yaml'

def specs():
    c=yaml.safe_load(CAMPAIGN.read_text())
    return [yaml.safe_load((ROOT/x['specification']).read_text()) for x in c['cases']]

def test_01_exact_signed_k_parsing(): assert exact_rational('-5/2')==QQ(-5)/2
def test_02_float_cs_rejected():
    with pytest.raises(ConventionError): exact_rational(-.5,'signed_k')
def test_03_unsigned_k_rejected():
    s=deepcopy(specs()[0]); s['absolute_k']='1/2'
    with pytest.raises(ConventionError): validate_spec_shape(s)
def test_04_BQ_normalization(): assert all(s['normalization']['B_of_Q']==1 for s in specs())
def test_05_beta_conversion(): assert beta_to_baryon(-1)==-3
def test_06_two_anchor_exact_solution():
    s=specs()[0]; M,_=solve_two_anchor_map(s['raw_charge_order'],s['classical_anchor'],s['instanton_anchor']); assert M==matrix(QQ,[[-QQ(5)/4,-QQ(7)/4],[QQ(1)/2,-QQ(1)/2]])
def test_07_rank_uniqueness_rejected():
    a={'raw_charges':{'x':1,'q':0},'target':{'B':3,'I':0}}
    with pytest.raises(ConventionError): solve_two_anchor_map(['x','q'],a,a)
def test_08_representation_aware_lookup():
    s=specs()[2]; raw=json.loads((ROOT/'generated/su3_nf5_k5o2_infinite/order_10/branching_comparison/raw_branching.json').read_text()); assert find_anchor(raw,s['instanton_anchor'])[1]['x_charge']==-4
def test_09_wrong_representation_rejected():
    s=deepcopy(specs()[2]); s['instanton_anchor']['child_representation']=[0,0,1,0]
    raw=json.loads((ROOT/'generated/su3_nf5_k5o2_infinite/order_10/branching_comparison/raw_branching.json').read_text())
    with pytest.raises(ConventionError): find_anchor(raw,s['instanton_anchor'])
def test_10_zero_mode_formula(): assert negative_cs_instanton_baryon(3,'-3/2')==-QQ(3)/2
def test_11_section_23_numeric(): assert negative_cs_instanton_baryon(3,'-1/2')==-QQ(5)/2
def test_12_operator_conjugation(): assert conjugate_operator({'representation':[1,0,0,0],'B':-2,'I':1})=={'representation':[0,0,0,1],'B':2,'I':-1}
def test_13_reverse_k_is_not_conjugation():
    a={'representation':[1,0,0,0],'B':-2,'I':1}; assert conjugate_operator(a)['I']==-1 and reverse_cs_orientation(a,'-1')['I']==1
def test_14_self_conjugate_20_beta_sign(): assert conjugate_labels([0,0,1,0,0])==[0,0,1,0,0] and beta_to_baryon(-1)==-3
def test_15_half_integral_5f_supported(): assert [s['charge_lattice']['B_step'] for s in specs()[:3]]==['1/2']*3
def test_16_integral_6f_supported(): assert [s['charge_lattice']['B_step'] for s in specs()[3:]]==['1']*3
def test_17_integer_I_all_cases(): assert all(validate_case(ROOT,ROOT/f'theories/branching/{s["theory_id"].replace("_infinite","")}_to_finite.yaml')['instanton_I']=='1' for s in specs())
def test_18_classical_I_zero(): assert all(s['classical_anchor']['target']['I']==0 for s in specs())
def test_19_current_neutral_rejected_and_classified():
    s=deepcopy(specs()[0]); s['instanton_anchor']['target']['B']=0
    assert classify_stored_map(0,1)=='legacy_current_neutral'
    with pytest.raises(ConventionError,match='legacy shifted basis'): validate_spec_shape(s)
def test_20_canonical_B_rendering():
    r=render_canonical_charge([1,0],-2,1); assert 'B=' in r and 'mic' not in r and 'widehat' not in r
def test_21_signed_title(): assert render_signed_title(3,'-3/2',5)==r'SU(3)_{-3/2}+5F'
def test_22_k5o2_conjugate_current(): assert specs()[2]['instanton_anchor']['raw_charges']=={'x':-4,'q':0}
def test_23_raw_branching_hashes_are_inputs():
    run_preflight(ROOT,CAMPAIGN); h=json.loads((ROOT/'generated/campaigns/su3_nf5_nf6_branching_preflight/input_hashes.json').read_text()); assert len([x for x in h if x.endswith('raw_branching.json')])==6
def test_24_deterministic_preflight():
    run_preflight(ROOT,CAMPAIGN); p=ROOT/'generated/campaigns/su3_nf5_nf6_branching_preflight'; before={x.name:x.read_bytes() for x in p.iterdir() if x.name!='command_log.txt'}; run_preflight(ROOT,CAMPAIGN); assert before=={x.name:x.read_bytes() for x in p.iterdir() if x.name!='command_log.txt'}
