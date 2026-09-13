# independent Python implementation of the reduced relaxation for cross-validation.
# Third argument 'nb': also build the restricted automaton R''_m (no close on an empty partial stack,
# no (close,close) popping a '?') and report its Perron root over the states reachable from the start.
import sys, numpy as np, scipy.sparse as sp
from collections import Counter
def canon(U,L):
    def lab(seq,m,nx):
        out=[]
        for x in seq:
            if x==0: out.append(0);continue
            if x not in m: m[x]=nx[0];nx[0]+=1
            out.append(m[x])
        return tuple(out)
    m={};nx=[1];a=(lab(U,m,nx),lab(L,m,nx))
    m={};nx=[1];b=(lab(L,m,nx),lab(U,m,nx))
    return min(a,b)
def strip(U,L):
    c=Counter([x for x in U+L if x])
    U=[x if (x and c[x]==2) else 0 for x in U];L=[x if (x and c[x]==2) else 0 for x in L]
    while U and U[0]==0: U.pop(0)
    while L and L[0]==0: L.pop(0)
    return U,L
def succ(state,k,TOT):
    U,L=state;out=[]
    for eu in (0,1):
        for el in (0,1):
            U2=list(U);L2=list(L);closed=[]
            under=(eu==0 and not U2) or (el==0 and not L2)
            if eu==0: closed.append(U2.pop() if U2 else 0)
            if el==0: closed.append(L2.pop() if L2 else 0)
            if len(closed)==2 and closed[0] and closed[0]==closed[1]: continue
            tag=1 if (under or (len(closed)==2 and (closed[0]==0 or closed[1]==0))) else 0
            cur=0
            if len(closed)==2:
                a,b=closed
                if a and b: U2=[a if x==b else x for x in U2];L2=[a if x==b else x for x in L2];cur=a
            elif closed: cur=closed[0]
            if eu or el:
                if not cur: cur=max(U2+L2+[0])+1
                if eu: U2.append(cur)
                if el: L2.append(cur)
            U2,L2=strip(U2,L2)
            if TOT>0:
                while len(U2)+len(L2)>TOT:
                    if len(U2)>=len(L2): U2.pop(0)
                    else: L2.pop(0)
            else:
                while len(U2)>k: U2.pop(0)
                while len(L2)>k: L2.pop(0)
            U2,L2=strip(U2,L2);U2,L2=strip(U2,L2)
            out.append((canon(tuple(U2),tuple(L2)),tag))
    return out
k=int(sys.argv[1]);TOT=int(sys.argv[2]);NB=len(sys.argv)>3 and sys.argv[3]=='nb'
start=canon((),());idx={start:0};states=[start];edges=[];tags=[];i=0
while i<len(states):
    for t,tag in succ(states[i],k,TOT):
        if t not in idx: idx[t]=len(states);states.append(t)
        edges.append((i,idx[t]));tags.append(tag)
    i+=1
n=len(states);A=sp.coo_matrix((np.ones(len(edges)),([e[0] for e in edges],[e[1] for e in edges])),shape=(n,n)).tocsr()
v=np.ones(n);ev=0
for it in range(3000):
    w=A@(A@v);m=w.max();nev=np.sqrt(m);v=w/m
    if abs(nev-ev)<1e-13: break
    ev=nev
cw=(A@(A@v)/v).max()
print(f"python k={k} tot={TOT} states={n} rho(T^2)~{ev*ev:.6f} CW upper {cw:.6f}")
if NB:
    keep=[e for e,t in zip(edges,tags) if not t]
    B=sp.coo_matrix((np.ones(len(keep)),([e[0] for e in keep],[e[1] for e in keep])),shape=(n,n)).tocsr()
    # reachable from start under B
    seen=np.zeros(n,bool);seen[0]=True;stack=[0]
    while stack:
        i=stack.pop()
        for j in B.indices[B.indptr[i]:B.indptr[i+1]]:
            if not seen[j]: seen[j]=True;stack.append(j)
    R=B[seen][:,seen];v=np.ones(R.shape[0]);ev=0
    for it in range(5000):
        w=R@(R@v);m=w.max();nev=np.sqrt(m);v=w/m
        if abs(nev-ev)<1e-13: break
        ev=nev
    r=R@(R@v)/v
    print(f"python restricted R''_{TOT}: reachable {seen.sum()} of {n}, rho(T''^2)~{ev*ev:.6f}, CW [{r.min():.6f},{r.max():.6f}]")
