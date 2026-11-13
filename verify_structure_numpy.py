#!/usr/bin/env python3
"""Full exact NumPy audit of every truth-table certificate.
Usage: python verify_structure_numpy.py [--root DIR] [--report FILE]
Requires Python 3.10+ and NumPy. Does not trust search logs or stored spectra.
"""
from __future__ import annotations
import argparse, hashlib, json, sys
from pathlib import Path
import numpy as np

def read_table(path:Path,m:int)->np.ndarray:
 if path.suffix=='.bin': return np.fromfile(path,dtype=np.uint8)
 text=path.read_text(encoding='ascii')
 text='\n'.join(s for s in text.splitlines() if not s.lstrip().startswith('#'))
 if m==1:
  s=''.join(text.split())
  if any(v not in '01' for v in s):raise ValueError('Boolean text contains non-binary symbols')
  return np.fromiter((ord(v)-48 for v in s),dtype=np.uint8)
 tokens=text.split()
 vals=[int(v) for v in tokens]
 if any(v<0 or v>=1<<m for v in vals):raise ValueError('Output label out of range')
 return np.array(vals,dtype=np.uint8)

def fwht(a:np.ndarray)->np.ndarray:
 v=np.array(a,dtype=np.int64,copy=True);h=1
 while h<len(v):
  z=v.reshape(-1,2*h);left=z[:,:h].copy();right=z[:,h:].copy()
  z[:,:h]=left+right;z[:,h:]=left-right;h*=2
 return v

def degree(f:np.ndarray,wt:np.ndarray)->int:
 v=np.array(f,dtype=np.uint8,copy=True);h=1
 while h<len(v):
  z=v.reshape(-1,2*h);z[:,h:]^=z[:,:h];h*=2
 return int(wt[v!=0].max()) if np.any(v) else 0

def expected_layers(entry:dict,l:int):
 n=entry['n'];N=1<<n;p=entry['profile'];kind=p['type']
 if kind=='field_trace':return None
 a=[0]*(n+1)
 if kind=='single_layer':a[p['j']]=N
 elif kind=='two_term':a[p['j1']]=2*N;a[p['j2']]=-N
 elif kind=='vector9_joint':
  if l in (1,5,7):a[5]=2*N;a[7]=-N
  elif l in (2,4,6):a[6]=N
  elif l==3:a[5]=-N;a[7]=2*N
  else:raise ValueError('Unknown component')
 else:raise ValueError('Unknown profile')
 return a

def verify(root:Path,entry:dict)->dict:
 path=root/entry['file'];n=entry['n'];m=entry['m'];N=1<<n
 raw=path.read_bytes();y=read_table(path,m)
 if len(y)!=N or np.any(y>=1<<m):raise ValueError('Wrong truth-table length or output range')
 wt=np.fromiter((x.bit_count() for x in range(N)),dtype=np.int16,count=N)
 cp={};checks=[]
 for l in range(1,1<<m):
  lut=np.array([(l&b).bit_count()&1 for b in range(1<<m)],dtype=np.uint8)
  f=lut[y];W=fwht(1-2*f.astype(np.int64))
  A=[int(W[wt==j].sum()) for j in range(n+1)]
  c=next(j for j,a in enumerate(A) if a)-1
  supp=W!=0;M=int(np.abs(W).max());res=int(wt[supp].min())-1
  expected=expected_layers(entry,l)
  cc=dict(c=c,M=M,NL=(N-M)//2,res=res,degree=degree(f,wt),balanced=bool(W[0]==0),
    layer_sums=A,parseval_ok=int(np.sum(W*W))==N*N,
    inverse_transform_ok=bool(np.array_equal(fwht(W),N*(1-2*f.astype(np.int64)))),
    four_divisibility_ok=bool(np.all(W%4==0)),
    required_layer_polynomial_ok=expected is None or A==expected)
  cp[str(l)]=cc
  checks.extend([cc['parseval_ok'],cc['inverse_transform_ok'],cc['four_divisibility_ok'],cc['required_layer_polynomial_ok']])
 pp=list(cp.values());actual=dict(c=min(p['c'] for p in pp),M=max(p['M'] for p in pp),NL=min(p['NL'] for p in pp),res=min(p['res'] for p in pp),degree=max(p['degree'] for p in pp),balanced=all(p['balanced'] for p in pp))
 comparison={k:actual[k]==v for k,v in entry['expected'].items()}
 counts=[[int(np.sum(y[wt==w]==b)) for b in range(1<<m)] for w in range(n+1)]
 uniform=all(int(np.sum(y==b))==N//(1<<m) for b in range(1<<m))
 return dict(file=entry['file'],n=n,m=m,sha256=hashlib.sha256(raw).hexdigest(),size_bytes=len(raw),expected=entry['expected'],actual=actual,parameter_checks=comparison,uniform_preimage_sizes=uniform,components=cp,preimage_layer_counts=counts,passed=all(checks) and all(comparison.values()) and uniform)

def main()->int:
 parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--root',type=Path,default=Path(__file__).resolve().parent);parser.add_argument('--report',type=Path)
 args=parser.parse_args();root=args.root.resolve();manifest=json.loads((root/'manifest.json').read_text());results=[]
 for e in manifest['entries']:
  try:r=verify(root,e)
  except Exception as exc:r=dict(file=e['file'],passed=False,error=f'{type(exc).__name__}: {exc}')
  results.append(r);print(('PASS' if r['passed'] else 'FAIL')+' '+e['file'],r.get('actual',r.get('error')))
 output=dict(verifier='Exact integer NumPy FWHT and Boolean Möbius inversion',all_passed=all(r['passed'] for r in results),total=len(results),passed=sum(r['passed'] for r in results),results=results)
 report=args.report or root/'verification/verification_numpy.json';report.parent.mkdir(parents=True,exist_ok=True);report.write_text(json.dumps(output,indent=2)+'\n')
 print(f"{output['passed']}/{output['total']} passed. Report: {report}")
 return 0 if output['all_passed'] else 1
if __name__=='__main__':raise SystemExit(main())
