from wcount import *
from sympy import Matrix, symbols, eye, cancel, factor
x=symbols('x')
def gf(W):
    from itertools import product
    def step(st):
        R,L=st; outs=[]
        for er,el in product([1,0],repeat=2):
            if not er and not R: continue
            if not el and not L: continue
            if len(R)+len(L)+(1 if er else -1)+(1 if el else -1)>W: continue
            R2=list(R);L2=list(L);touched=[]
            if not er: touched.append(R2.pop())
            if not el: touched.append(L2.pop())
            if len(touched)==2 and touched[0]==touched[1]:
                if not R2 and not L2: outs.append('DONE')
                continue
            if not touched: cid=max(list(R2)+list(L2)+[0])+1
            else:
                cid=min(touched); R2=[cid if c in touched else c for c in R2]; L2=[cid if c in touched else c for c in L2]
            if er: R2.append(cid)
            if el: L2.append(cid)
            outs.append(canon(tuple(R2),tuple(L2)))
        return outs
    init=((),()); states={init:0}; order=[init]; edges={}; i=0
    while i<len(order):
        s=order[i]; i+=1
        for t in step(s):
            if t!='DONE' and t not in states: states[t]=len(order); order.append(t)
            edges.setdefault(s,[]).append(t)
    N=len(order); M=Matrix(N,N,lambda i,j:0); d=Matrix(N,1,lambda i,j:0)
    for s,ts in edges.items():
        for t in ts:
            if t=='DONE': d[states[s]]+=1
            else: M[states[s],states[t]]+=1
    v=Matrix(1,N,lambda i,j:1 if j==0 else 0)
    return N,factor(cancel((v*(eye(N)-x*M).inv()*d)[0]*x))
for W in (4,6):
    N,G=gf(W); print("width<=%d states=%d GF:"%(W,N),G)
