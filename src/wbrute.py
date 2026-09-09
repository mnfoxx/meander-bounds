def matchings(pts):
    if not pts: return [[]]
    out=[]; a=pts[0]
    for j in range(1,len(pts),2):
        b=pts[j]
        for m1 in matchings(pts[1:j]):
            for m2 in matchings(pts[j+1:]): out.append([(a,b)]+m1+m2)
    return out
def hprof(m,N):
    h=[0]*(N+1)
    for a,b in m:
        for g in range(a,b): h[g]+=1
    return h
def single_cycle(up,lo,N):
    nu={};nl={}
    for a,b in up: nu[a]=b;nu[b]=a
    for a,b in lo: nl[a]=b;nl[b]=a
    cur=1;steps=0;side=0
    while True:
        cur=(nu if side==0 else nl)[cur]; steps+=1; side^=1
        if cur==1: return steps==N
for k in range(1,7):
    N=2*k; ms=matchings(list(range(1,N+1))); hp=[hprof(m,N) for m in ms]; byw={}
    for i,up in enumerate(ms):
        for j,lo in enumerate(ms):
            if single_cycle(up,lo,N):
                w=max(a+b for a,b in zip(hp[i],hp[j])); byw[w]=byw.get(w,0)+1
    print("k=%d width<=4,6,8:"%k,[sum(v for ww,v in byw.items() if ww<=W) for W in (4,6,8)])
