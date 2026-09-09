import sys, numpy as np, scipy.sparse as sp
from tilt import succ, canon
from collections import defaultdict
k=int(sys.argv[1]);TOT=int(sys.argv[2])
start=canon((),());idx={start:0};states=[start];rows=[];cols=[];i=0
while i<len(states):
    for t,dw in succ(states[i],k,TOT):
        if t not in idx: idx[t]=len(states);states.append(t)
        rows.append(i);cols.append(idx[t])
    i+=1
n=len(states);A=sp.coo_matrix((np.ones(len(rows)),(rows,cols)),shape=(n,n)).tocsr()
def perron(M):
    v=np.ones(M.shape[0]);ev=0
    for it in range(5000):
        w=M@(M@v);m=w.max();nev=np.sqrt(m);v=w/m
        if abs(nev-ev)<1e-13: break
        ev=nev
    # symmetrize to T eigenvector
    v=v+ (M@v)/ev
    return v/v.sum(),ev
b,ev=perron(A);a,_=perron(A.T.tocsr())
mass=a*b;mass/=mass.sum()
by=defaultdict(float);byq=defaultdict(float)
for s,mv in zip(states,mass):
    U,L=s;ne=len(U)+len(L);npair=len(set(x for x in U+L if x));nq=sum(1 for x in U+L if x==0)
    by[(ne,npair)]+=mv;byq[nq]+=mv
print("m",TOT,"rho^2",ev*ev)
print("mass by (#entries,#pairs):")
for key in sorted(by): 
    if by[key]>0.005: print("  ",key,f"{100*by[key]:.1f}%")
print("mass by #? entries:",{q:f"{100*v:.1f}%" for q,v in sorted(byq.items())})
