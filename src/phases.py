import sys
from sympy import Matrix, symbols, factor, Poly, eye
exec(open('tmD.py').read().split('import numpy')[0].replace('MAXH=int(sys.argv[1])','MAXH=int(sys.argv[1])').replace('print("states:",N)',''))
x=symbols('x')
def charpoly_rev(idx):
    m=len(idx); pos={s:i for i,s in enumerate(idx)}
    M=Matrix(m,m,lambda i,j:0)
    for s in idx:
        for t in edges.get(s,[]):
            if t!='DONE' and t in pos: M[pos[s],pos[t]]+=1
    return factor((eye(m)-x*M).det())
pre=[s for s in order if not s[2]]; post=[s for s in order if s[2]]
print("D=%d  states pre/post: %d/%d"%(MAXH-1,len(pre),len(post)))
# how many post states carry a two-ended component?
labelled=[s for s in post if any(c not in (0,'T') for c in s[0]+s[1])]
print("post states with two-ended paths:",len(labelled),"; without:",len(post)-len(labelled))
print("pre  charpoly:",charpoly_rev(pre))
print("post charpoly:",charpoly_rev(post))
