# Closed and open meanders of bounded height via automata; exact rational GFs.
import sys
from itertools import product
from fractions import Fraction
from sympy import Matrix, symbols, factor, eye, cancel, Poly, Rational
x=symbols('x')
def canon(R,L):
    m={};nxt=1
    def f(c):
        nonlocal nxt
        if c not in m: m[c]=nxt;nxt+=1
        return m[c]
    return (tuple(f(c) for c in R),tuple(f(c) for c in L))
def build_closed(H):
    # rootless automaton; accept = closing the two ends of the unique remaining path with stacks otherwise empty
    def step(st):
        R,L=st; outs=[]
        for er,el in product(['open','close'],repeat=2):
            if er=='close' and not R: continue
            if el=='close' and not L: continue
            if er=='open' and len(R)>=H: continue
            if el=='open' and len(L)>=H: continue
            R2=list(R);L2=list(L);touched=[]
            if er=='close': touched.append(R2.pop())
            if el=='close': touched.append(L2.pop())
            if len(touched)==2 and touched[0]==touched[1]:
                if not R2 and not L2: outs.append('DONE')
                continue
            if not touched: cid=max(list(R2)+list(L2)+[0])+1
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
            if t!='DONE' and t not in states: states[t]=len(order); order.append(t)
            edges.setdefault(s,[]).append(t)
    N=len(order); M=Matrix(N,N,lambda i,j:0); d=Matrix(N,1,lambda i,j:0)
    for s,ts in edges.items():
        for t in ts:
            if t=='DONE': d[states[s]]+=1
            else: M[states[s],states[t]]+=1
    v=Matrix(1,N,lambda i,j:1 if j==0 else 0)
    G=cancel((v*(eye(N)-x*M).inv()*d)[0]*x)   # x^(crossings): crossings = steps+1 (DONE step is a crossing)
    return N,G
for H in [2,3,4]:
    N,G=build_closed(H)
    ser=Poly(cancel(G).series(x,0,26).removeO(),x)
    coeffs=[ser.coeff_monomial(x**k) for k in range(26)]
    print("closed meanders, height<=%d: states=%d"%(H,N))
    print("  by crossings 2k, k=1..12:",[coeffs[2*k] for k in range(1,13)])
    print("  GF:",factor(G))
