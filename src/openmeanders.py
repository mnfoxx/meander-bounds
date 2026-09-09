import sys
from sympy import *
x=symbols('x')
for H in [2,3,4]:
    sys.argv=['x',str(H)]; ns={}
    exec(open('tmD.py').read().split('import numpy')[0].replace('print("states:",N)',''),ns)
    order=ns['order']; edges=ns['edges']; states=ns['states']
    pre=[s for s in order if not s[2]]; idx={s:i for i,s in enumerate(pre)}
    N=len(pre); M=Matrix(N,N,lambda i,j:0)
    for s in pre:
        for t in edges.get(s,[]):
            if t!='DONE' and t in idx: M[idx[s],idx[t]]+=1
    i0=idx[((0,),(),False)]; i1=idx[((),(0,),False)]
    Rm=(eye(N)-x*M).inv()
    # open meander with m crossings = folding of n=m+2 stamps, 1 top, n bottom: walk start->start (n even) or start->mirror (n odd), length n-2=m
    G=cancel(Rm[i0,i0]+Rm[i0,i1])
    ser=Poly(G.series(x,0,16).removeO(),x); co=[ser.coeff_monomial(x**k) for k in range(16)]
    print("open meanders height<=%d (m=0..15):"%H,co)
    print("   GF:",factor(G))
