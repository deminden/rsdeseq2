#!/usr/bin/env python3
"""Independent 90-decimal KKT roots for frozen original mixture likelihoods.

No candidate executable or candidate numerical output is read. Enumerating active
faces and certifying the full concave likelihood gradient establishes the global
maximum, independently of either the native barrier solver or original SQP.
"""
import argparse, itertools, json, math, struct
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import mpmath as mp

def exact(value):
    if isinstance(value,str) and value.startswith('bits:'):
        value=struct.unpack('>d',bytes.fromhex(value[5:]))[0]
    return mp.mpf(value)

def solve(c):
    n,k=c['observations'],c['components']
    ll=list(map(exact,c['log_likelihoods']))
    a=list(map(exact,c['observation_weights'])) if c['observation_weights'] else [mp.mpf(1)]*n
    prior=[exact(v)-1 for v in c['dirichlet']];mass=sum(a)+sum(prior)
    maxima=[max(ll[i*k:(i+1)*k])for i in range(n)]
    L=[[mp.exp(ll[i*k+j]-maxima[i]) for j in range(k)]for i in range(n)]
    initial=list(map(exact,c['initial_weights']))
    mandatory={j for j,p in enumerate(prior)if p>0}
    for size in range(len(mandatory),k+1):
        for active in itertools.combinations(range(k),size):
            if not mandatory.issubset(active):continue
            def gradient(*p):
                marginal=[sum(row[j]*v for j,v in zip(active,p))for row in L]
                return tuple(sum(a[i]*L[i][j]/marginal[i]for i in range(n))+prior[j]/p[t]-mass for t,j in enumerate(active))
            def jacobian(*p):
                marginal=[sum(row[j]*v for j,v in zip(active,p))for row in L]
                return mp.matrix([[-sum(a[i]*L[i][j]*L[i][l]/marginal[i]**2 for i in range(n))-(prior[j]/p[t]**2 if j==l else 0)for l in active]for t,j in enumerate(active)])
            start=[max(initial[j],mp.mpf('1e-3'))for j in active];total=sum(start);start=[v/total for v in start]
            try:p=list(mp.findroot(gradient,start,J=jacobian,tol=mp.mpf('1e-75'),maxsteps=100))
            except (ValueError,ZeroDivisionError):continue
            if any(v<=0 for v in p):continue
            if abs(sum(p)-1)>mp.mpf('1e-65'):continue
            full=[mp.mpf(0)]*k
            for j,v in zip(active,p):full[j]=v
            marginal=[sum(row[j]*v for j,v in enumerate(full))for row in L]
            g=[sum(a[i]*L[i][j]/marginal[i]for i in range(n))+(prior[j]/full[j]if prior[j]>0 else 0)for j in range(k)]
            if any(g[j]>mass+mp.mpf('1e-60')for j in range(k)if j not in active):continue
            likelihood=sum(a[i]*(maxima[i]+mp.log(marginal[i]))for i in range(n))
            penalty=sum(prior[j]*mp.log(full[j])for j in active if prior[j]>0)
            return dict(name=c['name'],weights=[mp.nstr(v,82)for v in full],log_likelihood=mp.nstr(likelihood,82),penalized_log_likelihood=mp.nstr(likelihood+penalty,82),active=list(active),mass=mp.nstr(mass,82),gradient=[mp.nstr(v,82)for v in g])
    raise RuntimeError('No independently certified active face for '+c['name'])

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('reference',type=Path);parser.add_argument('output',type=Path)
    args=parser.parse_args()
    if args.output.exists():raise SystemExit('Use a fresh output file')
    mp.mp.dps=90
    source=json.loads(args.reference.read_text())
    cases=[solve(c)for c in source['em_fit_cases']if c['optimizer']=='mixIP']
    args.output.write_text(json.dumps(dict(source='90-decimal active-face KKT roots of frozen original likelihoods',mpmath=mp.__version__,cases=cases),indent=2)+'\n')

if __name__=='__main__':main()
