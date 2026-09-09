import sys
sys.setrecursionlimit(10000)
from irred import matchings, is_meander
from collections import Counter
def splits(up,lo):
    n=len(up);res=[]
    for c in range(1,n):
        for side,m in (('U',up),('L',lo)):
            if m[c-1]==c:
                sp=[i for i in range(c) if m[i]>=c]
                other=lo if side=='U' else up
                if len(sp)==2 and not any(other[i]>=c for i in range(c)):
                    res.append((c,side))
    return res
def undo(up,lo,c,side):
    up=up[:];lo=lo[:];m=up if side=='U' else lo
    sp=sorted(i for i in range(c) if m[i]>=c);p=sp[0];q=m[p]
    assert sp[1]==c-1 and m[c-1]==c and p<c-1<c<q
    m[p]=c-1;m[c-1]=p;m[c]=q;m[q]=c
    return up,lo
def factor(up,lo):
    # undo all split positions (they are determined); check independence of order
    S=splits(up,lo)
    for c,side in S: up,lo=undo(up,lo,c,side)
    assert splits(up,lo)==[] or True
    n=len(up);seen=[False]*n;bl=[]
    for s in range(n):
        if seen[s]:continue
        comp=[];i=s
        while not seen[i]:
            seen[i]=True;comp.append(i);j=up[i];seen[j]=True;comp.append(j);i=lo[j]
        comp.sort();assert comp==list(range(comp[0],comp[-1]+1));bl.append(comp)
    keys=[]
    for comp in bl:
        o=comp[0];keys.append((tuple(up[i]-o for i in comp),tuple(lo[i]-o for i in comp)))
    return keys,S
K=int(sys.argv[1]);M={};P={}
for k in range(1,K+1):
    ms=matchings(2*k);tot=0;pr=0;bad=0;comp=Counter()
    for up in ms:
        for lo in ms:
            if is_meander(up,lo):
                tot+=1
                keys,S=factor(up,lo)
                if not S: pr+=1
                for key in keys:
                    if splits(list(key[0]),list(key[1])): bad+=1
                comp[(tuple(len(x[0])//2 for x in keys),tuple(s for _,s in S))]+=1
    M[k]=tot;P[k]=pr
    print(k,tot,pr,"nonprime blocks:",bad)
    if k<=5:
        for key,v in sorted(comp.items()): print("   ",key,v)
# test M = 1/(1-2P)? sequences with side choice: each join has 2 sides -> M = 1/(1 - P) with P counting... check numerically
import fractions
def seq(Pw):
    # coefficients of 1/(1-Pw(x))
    Mc=[1]
    for n in range(1,K+1):
        Mc.append(sum(Pw.get(j,0)*Mc[n-j] for j in range(1,n+1)))
    return Mc
print("M:",[M[k] for k in range(1,K+1)])
for name,Pw in (("1/(1-P)",P),("1/(1-2P)",{k:2*v for k,v in P.items()})):
    print(name,seq(Pw)[1:])
