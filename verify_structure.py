#!/usr/bin/env python3
"""Independent, standard-library-only verification of the appendix.

Unlike verify_structure.py, this verifier uses scalar Python integer
butterflies and computes layer polynomials independently from preimage counts
and Krawtchouk coefficients. No NumPy or generation module is imported.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
from pathlib import Path

def read_outputs(path:Path,m:int)->list[int]:
 if path.suffix=='.bin':return list(path.read_bytes())
 lines=[line for line in path.read_text(encoding='ascii').splitlines() if not line.lstrip().startswith('#')]
 text='\n'.join(lines)
 if m==1:
  text=''.join(text.split())
  if any(c not in '01' for c in text):raise ValueError('Invalid Boolean character')
  return [int(c) for c in text]
 return [int(c) for c in text.split()]

def hadamard(values:list[int])->list[int]:
 out=values.copy();N=len(out);stride=1
 while stride<N:
  for start in range(0,N,2*stride):
   for offset in range(stride):
    i=start+offset;j=i+stride;left=out[i];right=out[j]
    out[i]=left+right;out[j]=left-right
  stride*=2
 return out

def algebraic_degree(outputs:list[int],wt:list[int])->int:
 # Möbius transformation of packed output words checks all coordinates together.
 a=outputs.copy();N=len(a);bit=1
 while bit<N:
  for base in range(0,N,2*bit):
   for j in range(bit):a[base+bit+j]^=a[base+j]
  bit*=2
 return max((wt[i] for i,v in enumerate(a) if v),default=0)

def krawtchouk(n:int,j:int,w:int)->int:
 # Degree j, evaluated at input weight w.
 return sum((-1 if r&1 else 1)*math.comb(w,r)*math.comb(n-w,j-r)
            for r in range(max(0,j-(n-w)),min(w,j)+1))

def expected_polynomial(entry:dict,l:int):
 n=entry['n'];N=1<<n;p=entry['profile'];a=[0]*(n+1)
 if p['type']=='field_trace':return None
 if p['type']=='single_layer':a[p['j']]=N
 elif p['type']=='two_term':a[p['j1']]=2*N;a[p['j2']]=-N
 elif p['type']=='vector9_joint':
  if l==3:a[5]=-N;a[7]=2*N
  elif l&1:a[5]=2*N;a[7]=-N
  else:a[6]=N
 else:raise ValueError('Unknown profile')
 return a

def multiply(a:int,b:int,n:int,poly:int)->int:
 r=0
 while b:
  if b&1:r^=a
  b>>=1;a<<=1
  if a>>n:a^=poly
 return r

def power(a:int,e:int,n:int,poly:int)->int:
 r=1
 while e:
  if e&1:r=multiply(r,a,n,poly)
  a=multiply(a,a,n,poly);e//=2
 return r

def trace(a:int,n:int,poly:int)->int:
 r=0
 for _ in range(n):r^=a;a=multiply(a,a,n,poly)
 if r not in (0,1):raise ValueError('Trace is not binary')
 return r

def field_check(y:list[int],info:dict)->bool:
 n=info['n'];p=info['irreducible_polynomial_integer'];basis=info['ordered_input_basis'];gamma=info['trace_coefficients'];e=info['exponent']
 # Confirm that the field polynomial supports Frobenius x^(2^n)=x and
 # that the chosen input basis and output linear combinations are independent.
 if any(power(x,1<<n,n,p)!=x for x in range(1<<n)):return False
 def span(vectors):
  s={0}
  for v in vectors:s|={a^v for a in list(s)}
  return s
 if len(span(basis))!=1<<n or len(span(gamma))!=1<<len(gamma):return False
 for x,actual in enumerate(y):
  X=0
  for i,b in enumerate(basis):
   if x>>i&1:X^=b
  z=power(X,e,n,p);expected=0
  for i,g in enumerate(gamma):expected|=trace(multiply(g,z,n,p),n,p)<<i
  if actual!=expected:return False
 return True

def verify(root:Path,e:dict,field_info:dict)->dict:
 path=root/e['file'];n=e['n'];m=e['m'];N=1<<n;q=1<<m;y=read_outputs(path,m)
 if len(y)!=N or any(v<0 or v>=q for v in y):raise ValueError('Invalid truth table')
 wt=[x.bit_count() for x in range(N)];counts=[[0]*q for _ in range(n+1)]
 for x,v in enumerate(y):counts[wt[x]][v]+=1
 components={};global_M=0;global_c=n;global_res=n;all_checks=True
 for l in range(1,q):
  sign=[1-2*((l&v).bit_count()&1) for v in y];W=hadamard(sign)
  M=max(map(abs,W));A=[0]*(n+1)
  res=n
  for u,z in enumerate(W):
   A[wt[u]]+=z
   if z:res=min(res,wt[u]-1)
  S=[sum((1-2*((l&b).bit_count()&1))*counts[w][b] for b in range(q)) for w in range(n+1)]
  AK=[sum(S[w]*krawtchouk(n,j,w) for w in range(n+1)) for j in range(n+1)]
  if A!=AK:raise ValueError('Walsh layers differ from independent Krawtchouk moments')
  c=next(j for j,a in enumerate(AK) if a)-1
  expected=expected_polynomial(e,l)
  ok=(W[0]==0 and sum(z*z for z in W)==N*N and (expected is None or A==expected))
  all_checks=all_checks and ok
  components[str(l)]=dict(c=c,M=M,res=res,layer_sums=AK,balanced=W[0]==0,checks_passed=ok)
  global_M=max(global_M,M);global_c=min(global_c,c);global_res=min(global_res,res)
 actual=dict(c=global_c,M=global_M,NL=(N-global_M)//2,res=global_res,degree=algebraic_degree(y,wt),balanced=all(p['balanced'] for p in components.values()))
 field_ok=None
 if path.name in field_info:field_ok=field_check(y,field_info[path.name]);all_checks=all_checks and field_ok
 return dict(file=e['file'],sha256=hashlib.sha256(path.read_bytes()).hexdigest(),actual=actual,expected=e['expected'],components=components,krawtchouk_crosscheck_passed=True,field_realization_passed=field_ok,passed=all_checks and actual==e['expected'])

def main()->int:
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,default=Path(__file__).resolve().parent);p.add_argument('--report',type=Path);args=p.parse_args();root=args.root.resolve()
 manifest=json.loads((root/'manifest.json').read_text());infos=json.loads((root/'verification/field_constructions.json').read_text());field_info={x['file']:x for x in infos};results=[]
 for e in manifest['entries']:
  try:r=verify(root,e,field_info)
  except Exception as exc:r=dict(file=e['file'],passed=False,error=f'{type(exc).__name__}: {exc}')
  results.append(r);print(('PASS' if r['passed'] else 'FAIL')+' '+e['file'],r.get('actual',r.get('error')),flush=True)
 report=args.report or root/'verification/verification_independent.json'
 out=dict(verifier='Standard-library scalar integer FWHT, packed-output Möbius inversion and independent Krawtchouk moments',all_passed=all(x['passed'] for x in results),total=len(results),passed=sum(x['passed'] for x in results),results=results)
 report.parent.mkdir(parents=True,exist_ok=True);report.write_text(json.dumps(out,indent=2)+'\n');print(f"{out['passed']}/{out['total']} passed. Report: {report}")
 return 0 if out['all_passed'] else 1
if __name__=='__main__':raise SystemExit(main())
