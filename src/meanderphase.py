# "all two-ended paths" automaton: pairs of bounded stacks, labels = path ids, no cycles, no absents.
import sys
from itertools import product
from sympy import Matrix, symbols, factor, eye
MAXH=int(sys.argv[1])
def canon(R,L):
    m={};nxt=1
    def f(c):
        nonlocal nxt
        if c not in m: m[c]=nxt;nxt+=1
        return m[c]
    return (tuple(f(c) for c in R),tuple(f(c) for c in L))
def step(st):
    R,L=st; outs=[]
    for er,el in product(['open','close'],repeat=2):
        if er=='close' and not R: continue
        if el=='close' and not L: continue
        if er=='open' and len(R)>=MAXH: continue
        if el=='open' and len(L)>=MAXH: continue
        R2=list(R);L2=list(L);touched=[]
        if er=='close': touched.append(R2.pop())
        if el=='close': touched.append(L2.pop())
        if len(touched)==2 and touched[0]==touched[1]: continue  # cycle
        if not touched:
            cid=max(list(R2)+list(L2)+[0])+1
        else:
            cid=min(touched); R2=[cid if c in touched else c for c in R2]; L2=[cid if c in touched else c for c in L2]
        if er=='open': R2.append(cid)
        if el=='open': L2.append(cid)
        outs.append(canon(tuple(R2),tuple(L2)))
    return outs
init=((),()); states={init:0}; order=[init]; edges={}; i=0
while i<len(order):
    s=order[i]; i+=1
    for t in step(s):
        if t not in states: states[t]=len(order); order.append(t)
        edges.setdefault(s,[]).append(t)
N=len(order); x=symbols('x')
M=Matrix(N,N,lambda i,j:0)
for s,ts in edges.items():
    for t in ts: M[states[s],states[t]]+=1
print("D=%d meander-phase states: %d"%(MAXH-1,N))
print("charpoly:",factor((eye(N)-x*M).det()))
