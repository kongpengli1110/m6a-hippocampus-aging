"""Offline engineering tests of the actual aggregation helper, extracted by AST.
No biological analysis; exercises sparse grouping, repeated sample labels,
zero cells, noncontiguous groups, explicit zeros and noncanonical CSR inputs.
"""

# Release-only safety/path adapter; see provenance/adaptation_diffs.
from pathlib import Path as _ReleasePath
import sys as _release_sys
_release_home = next(p for p in _ReleasePath(__file__).resolve().parents if (p / 'release_support.py').is_file())
_release_sys.path.insert(0, str(_release_home))
from release_support import gate, release_root, resolve_source
gate(__file__, __name__)
from pathlib import Path
import ast,json,hashlib,datetime
import numpy as np,scipy.sparse as sp
root=Path(__file__).resolve().parents[1]
p=root/'scripts/03_aggregate_and_qa.py';tree=ast.parse(p.read_text())
f=[n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='accumulate_chunk'];assert len(f)==1
scope={'np':np,'sp':sp};exec(compile(ast.Module(body=f,type_ignores=[]),str(p),'exec'),scope)
accumulate=scope['accumulate_chunk'];rng=np.random.default_rng(20261012)
tested=0
for cells,genes,groups in [(1,1,1),(10,7,15),(513,109,37),(1025,199,80)]:
 x=rng.poisson(.6,size=(cells,genes)).astype(np.int64);x[rng.random(x.shape)<.7]=0
 labels=rng.integers(0,groups,size=cells,dtype=np.int32)
 expect=np.zeros((groups,genes),dtype=np.int64)
 for i in range(cells):expect[labels[i]]+=x[i]
 for step in [1,7,512,1024]:
  out=np.zeros_like(expect)
  for lo in range(0,cells,step):accumulate(sp.csr_matrix(x[lo:lo+step]),labels[lo:lo+step],out)
  assert np.array_equal(out,expect);tested+=1
# CSR with repeated coordinates, unsorted columns and explicit zero.
b=sp.csr_matrix((np.array([4,0,3,2,7],dtype=np.int64),np.array([2,0,2,1,1]),np.array([0,3,5])),shape=(2,4))
o=np.zeros((5,4),dtype=np.int64);accumulate(b,np.array([4,4]),o)
assert np.array_equal(o[4],np.array([0,9,7,0]));assert o[:4].sum()==0;tested+=1
report={'status':'PASS','tests':tested,'seed':20261012,'script_sha256_at_test':hashlib.sha256(p.read_bytes()).hexdigest(),'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
(root/'logs/synthetic_unit_tests.json').write_text(json.dumps(report,indent=2));print(report)
