import sys
from functools import lru_cache
sys.setrecursionlimit(10000)

def matchings(n):  # non-crossing matchings on points 0..n-1 as partner lists
    res=[]
    def rec(pts,cur):
        if not pts: res.append(cur[:]);return
        a=pts[0]
        for i in range(1,len(pts),2):
            b=pts[i]
            cur[a]=b;cur[b]=a
            rec(pts[1:i],cur);  # inside must be matched internally -> handled by splitting
    # simpler: standard recursion producing full partner arrays
    out=[]
    def gen(pts):
        if not pts: yield {}; return
        a=pts[0]
        for i in range(1,len(pts),2):
            b=pts[i]
            for m1 in gen(pts[1:i]):
                for m2 in gen(pts[i+1:]):
                    m=dict(m1);m.update(m2);m[a]=b;m[b]=a
                    yield m
    for m in gen(list(range(n))):
        out.append([m[i] for i in range(n)])
    return out

def is_meander(up,lo):
    n=len(up);seen=[False]*n;i=0;cnt=0
    while not seen[i]:
        seen[i]=True;cnt+=1
        j=up[i];seen[j]=True;cnt+=1;i=lo[j]
    return cnt==n

def cuts(up,lo):
    n=len(up);res=[]
    for c in range(1,n):  # cut between c-1 and c
        u=sum(1 for i in range(c) if up[i]>=c)
        l=sum(1 for i in range(c) if lo[i]>=c)
        if (u==2 and l==0) or (u==0 and l==2): res.append((c,'U' if u==2 else 'L'))
    return res

K=int(sys.argv[1])
for k in range(1,K+1):
    ms=matchings(2*k)
    tot=0;irr=0;bycuts={}
    for up in ms:
        for lo in ms:
            if is_meander(up,lo):
                tot+=1;c=cuts(up,lo)
                if not c: irr+=1
                bycuts[len(c)]=bycuts.get(len(c),0)+1
    print(k,tot,irr,sorted(bycuts.items()),flush=True)

def undo_all(up,lo):
    up=up[:];lo=lo[:]
    for c,side in cuts(up,lo):
        m=up if side=='U' else lo
        sp=[i for i in range(c) if m[i]>=c]
        u1,u2=sorted(sp)  # u1<u2<c ; w2=m[u2]<w1=m[u1]
        w1,w2=m[u1],m[u2]
        assert u1<u2<w2<w1
        m[u1]=u2;m[u2]=u1;m[w1]=w2;m[w2]=w1
    return up,lo

def blocks(up,lo):
    up2,lo2=undo_all(up,lo);n=len(up)
    # components
    seen=[False]*n;bl=[]
    for s in range(n):
        if seen[s]:continue
        comp=[];i=s
        while not seen[i]:
            seen[i]=True;comp.append(i);j=up2[i];seen[j]=True;comp.append(j);i=lo2[j]
        comp.sort();bl.append(comp)
    # each comp must be contiguous
    for comp in bl: assert comp==list(range(comp[0],comp[-1]+1)),comp
    # canonical block: relabel
    keys=[]
    for comp in bl:
        o=comp[0];u=tuple(up2[i]-o for i in comp);l=tuple(lo2[i]-o for i in comp)
        keys.append((u,l))
    return keys

from collections import Counter
if len(sys.argv)>2:
    for k in range(2,K+1):
        ms=matchings(2*k);cnt=Counter();irrset=set()
        for up in ms:
            for lo in ms:
                if is_meander(up,lo):
                    b=blocks(up,lo)
                    for key in b:
                        # check block irreducible
                        assert not cuts(list(key[0]),list(key[1])), key
                    cnt[tuple(len(key[0])//2 for key in b)]+=1
        # join-count per composition normalised by product of I's
        I={1:1,5:4,6:52,7:562}
        for comp,v in sorted(cnt.items()):
            prod=1
            for a in comp: prod*=I.get(a,0)
            print(k,comp,v, v/prod if prod else None)
