import sys
from sympy import *
exec(open('tmD.py').read().split('import numpy')[0].replace('print("states:",N)',''))
x=symbols('x')
post=[s for s in order if s[2]]
# recurrent post states: those on a cycle (reach themselves)
import collections
adj={s:[t for t in edges.get(s,[]) if t!='DONE'] for s in order}
def reach(s):
    seen=set(); st=list(adj[s])
    while st:
        u=st.pop()
        if u in seen: continue
        seen.add(u); st.extend(adj[u])
    return seen
rec=[s for s in post if s in reach(s)]
print("D=%d recurrent post-q states: %d"%(MAXH-1,len(rec)))
for s in rec: print("  R=%s L=%s"%(s[0],s[1]))
# tail GF: from the canonical state T=((0,),('T',),True), number of ways to DONE
T=((0,),('T',),True); idx={s:i for i,s in enumerate(order)}
N=len(order); M=Matrix(N,N,lambda i,j:0); d=Matrix(N,1,lambda i,j:0)
for s,ts in edges.items():
    for t in ts:
        if t=='DONE': d[idx[s]]+=1
        else: M[idx[s],idx[t]]+=1
v=Matrix(1,N,lambda i,j:1 if j==idx[T] else 0)
G=cancel((v*(eye(N)-x*M).inv()*d)[0])
print("tail GF from T:",factor(G))
