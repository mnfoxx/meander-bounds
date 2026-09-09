from sympy import Matrix, Rational, symbols, factor, nsimplify
from itertools import product
def guess_precursive(seq, maxorder=4, maxdeg=3):
    # find sum_{i=0..r} P_i(n) a(n-i) = 0 with deg P_i <= d, minimal (r,d) with nullspace
    N=len(seq)
    for r in range(1,maxorder+1):
        for d in range(0,maxdeg+1):
            unknowns=(r+1)*(d+1); eqs=N-r
            if eqs<unknowns+2: continue
            rows=[]
            for n in range(r,N):
                row=[]
                for i in range(r+1):
                    for j in range(d+1):
                        row.append(Rational(n)**j*seq[n-i])
                rows.append(row)
            M=Matrix(rows); ns=M.nullspace()
            if ns:
                v=ns[0]; v=v/max(abs(x) for x in v if x!=0)
                return r,d,list(v)
    return None
def poly_fit(seq):
    # check finite differences
    s=list(seq)
    for k in range(8):
        if len(set(s))==1: return k
        s=[s[i+1]-s[i] for i in range(len(s)-1)]
        if len(s)<3: return None
    return None
import sys
data={
 'depth1':[2,2,8,16,34,60,112,190,336,560,962,1588,2680,4398,7336,11992],
 'depth2':[0,0,0,6,30,94,289,732,1935,4547,11226,25468,60557,134802,313133,689300],
 'hair2':[1,5,6,11,14,19,22,31,34,45,50,57,64,79,84],
 'hair3':[1,12,26,49,99,158,256,371,561,750,1100,1401,1901,2340],
 'bottom4':[2,4,10,16,42,108,304,836,2454,7046,21276,62956,194008,586970],
}
n=symbols('n')
for name,s in data.items():
    print(name, 'poly degree(all):',poly_fit(s),' even:',poly_fit(s[0::2]),' odd:',poly_fit(s[1::2]))
    g=guess_precursive(s)
    print('  P-recursive guess:',g)
