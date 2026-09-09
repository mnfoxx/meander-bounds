import sys
from functools import lru_cache
def matchings(pts):
    # all non-crossing perfect matchings on list pts (sorted), as list of (a,b)
    if not pts: return [[]]
    out=[]; a=pts[0]
    for j in range(1,len(pts),2):
        b=pts[j]; inner=pts[1:j]; outer=pts[j+1:]
        for m1 in matchings(inner):
            for m2 in matchings(outer): out.append([(a,b)]+m1+m2)
    return out
def height(m,N):
    h=[0]*(N+1)
    for a,b in m:
        for g in range(a,b): h[g]+=1   # gap g between point g and g+1 (1-based points)
    return max(h)
def single_cycle(up,lo,N):
    nxt={}
    for a,b in up: nxt.setdefault(a,[]).append(b); nxt.setdefault(b,[]).append(a)
    for a,b in lo: nxt.setdefault(a,[]).append(b); nxt.setdefault(b,[]).append(a)
    # walk
    seen={1}; prev=None; cur=1; steps=0
    # alternate: start along upper arc
    side=0; 
    while True:
        arcs=up if side==0 else lo
        for a,b in arcs:
            if a==cur: nx=b; break
            if b==cur: nx=a; break
        cur=nx; steps+=1; side^=1
        if cur==1: return steps==N
        seen.add(cur)
for k in range(1,8):
    N=2*k; pts=list(range(1,N+1)); ms=matchings(pts)
    byh={}
    for up in ms:
        hu=height(up,N)
        for lo in ms:
            if single_cycle(up,lo,N):
                h=max(hu,height(lo,N)); byh[h]=byh.get(h,0)+1
    cum=[sum(v for hh,v in byh.items() if hh<=H) for H in (2,3,4)]
    print("k=%d total=%d  height<=2,3,4: %s"%(k,sum(byh.values()),cum))
