# (1) are all non-crossing cut-line configurations reachable? (2) shortest-path length vs width
from collections import deque
def canon(U,L):
    m={};nx=1;out=[]
    for x in U+L:
        if x not in m: m[x]=nx;nx+=1
    return (tuple(m[x] for x in U),tuple(m[x] for x in L))
def succ(state,W):
    U,L=state;out=[]
    for eu in (0,1):
        for el in (0,1):
            U2=list(U);L2=list(L);closed=[]
            if eu==0:
                if not U2: continue
                closed.append(U2.pop())
            if el==0:
                if not L2: continue
                closed.append(L2.pop())
            if len(U2)+len(L2)+eu+el>W: continue
            if len(closed)==2 and closed[0]==closed[1]: continue
            if len(closed)==2:
                a,b=closed;U2=[a if x==b else x for x in U2];L2=[a if x==b else x for x in L2];cur=a
            elif closed: cur=closed[0]
            else: cur=max(U2+L2+[0])+1
            if eu: U2.append(cur)
            if el: L2.append(cur)
            out.append(canon(tuple(U2),tuple(L2)))
    return out
W=10
dist={((),()):0};q=deque([((),())])
while q:
    s=q.popleft()
    for t in succ(s,W):
        if t not in dist: dist[t]=dist[s]+1;q.append(t)
import collections
byw=collections.defaultdict(list)
for s,d in dist.items(): byw[len(s[0])+len(s[1])].append(d)
for w in sorted(byw): print("width",w,"states",len(byw[w]),"min dist",min(byw[w]),"max dist",max(byw[w]),"max dist/width",max(byw[w])/w if w else 0)
# all non-crossing configurations of width w: matchings on w points on a line x split r
def matchings(n):
    if n==0: yield []; return
    for i in range(1,n,2):
        for m1 in matchings(i-1):
            for m2 in matchings(n-i-1):
                yield [(0,i)]+[(a+1,b+1) for a,b in m1]+[(a+i+1,b+i+1) for a,b in m2]
reach=set(dist)
for w in range(2,W+1,2):
    tot=0;ok=0
    for m in matchings(w):
        for r in range(0,w+1):
            # positions 0..r-1 = U outer->inner, r..w-1 = L inner->outer
            lab={}
            for k,(a,b) in enumerate(m): lab[a]=k+1;lab[b]=k+1
            U=tuple(lab[i] for i in range(r)); L=tuple(lab[i] for i in range(w-1,r-1,-1))
            tot+=1
            if canon(U,L) in reach: ok+=1
    print("width",w,"configs",tot,"reachable",ok)
