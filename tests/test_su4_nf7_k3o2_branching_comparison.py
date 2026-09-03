import json
from pathlib import Path
ROOT=Path(__file__).parents[1]
def test_complete_comparison_and_low_degrees():
 p=json.loads((ROOT/'generated/su4_nf7_k3o2_infinite/order_10/branching_comparison/finite_uv_comparison.json').read_text())
 assert [x['degree'] for x in p['degree_summary']]==[2,4,6,8,10]
 assert sum(p['summary'].values())==len(p['finite_terms'])
 by={(x['degree'],tuple(x['labels']),x['B']):x for x in p['finite_terms']}
 assert by[(2,(1,0,0,0,0,1),0)]['status']=='exact-match'
 assert by[(2,(0,0,0,0,0,0),0)]['uv_signed_multiplicity']==2
 assert by[(4,(0,0,1,0,0,0),-4)]['status']=='exact-match'
 assert by[(4,(0,0,0,0,0,0),0)]['uv_signed_multiplicity']==-2
