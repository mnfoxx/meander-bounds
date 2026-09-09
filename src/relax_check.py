# Mechanical check of Lemma L1 (soundness of the label-forgetting relaxation): explores all pairs
# (relaxed state r, exact state s) with r covering s, reachable from (empty, empty) with exact width <= W,
# and verifies that every exact transition is simulated by the relaxed automaton with covering preserved.
# Usage: relax_check.py W k [TOT]   (per-side memory k, or total memory TOT)
import sys

# relaxed successor keyed by event: replicate order: rsucc iterates eu in (0,1), el in (0,1) but skips illegal; rebuild with events
def rs(state,k):
    out={}
    U,L=state
    for eu in (0,1):
        for el in (0,1):
            U2=list(U);L2=list(L);closed=[];ok=True
            for side,e,S in ((0,eu,U2),(1,el,L2)):
                if e==0:
                    if S:
                        x=S.pop(); closed.append(x if x else None)
                    else: closed.append(None)
            if len(closed)==2 and closed[0] is not None and closed[0]==closed[1]: continue
            cur=None
            if len(closed)==2:
                a,b=closed
                if a is not None and b is not None:
                    U2=[a if x==b else x for x in U2];L2=[a if x==b else x for x in L2];cur=a
            elif len(closed)==1: cur=closed[0]
            if eu or el:
                if cur is None: cur=max(U2+L2+[0])+1
                if eu: U2.append(cur)
                if el: L2.append(cur)
            def strip(U2,L2):
                from collections import Counter
                c=Counter(U2+L2)
                U2=[x if (x and c[x]==2) else 0 for x in U2];L2=[x if (x and c[x]==2) else 0 for x in L2]
                while U2 and U2[0]==0: U2.pop(0)
                while L2 and L2[0]==0: L2.pop(0)
                return U2,L2
            U2,L2=strip(U2,L2)
            if TOT>0:
                while len(U2)+len(L2)>TOT:
                    if len(U2)>=len(L2): U2.pop(0)
                    else: L2.pop(0)
            else:
                for S in (U2,L2):
                    while len(S)>k: S.pop(0)
            U2,L2=strip(U2,L2);U2,L2=strip(U2,L2)
            out[(eu,el)]=(tuple(U2),tuple(L2))
    return out
def es(state,W):
    out={}
    U,L=state
    for eu in (0,1):
        for el in (0,1):
            U2=list(U);L2=list(L);closed=[]
            if eu==0:
                if not U2: continue
                closed.append(U2.pop())
            if el==0:
                if not L2: continue
                closed.append(L2.pop())
            if len(U2)+len(L2)+(eu+el)>W: continue
            if len(closed)==2 and closed[0]==closed[1]: continue
            if len(closed)==2:
                a,b=closed;U2=[a if x==b else x for x in U2];L2=[a if x==b else x for x in L2];cur=a
            elif closed: cur=closed[0]
            else: cur=max(U2+L2+[0])+1
            if eu: U2.append(cur)
            if el: L2.append(cur)
            out[(eu,el)]=(tuple(U2),tuple(L2))
    return out
def covers(r,s):
    (rU,rL),(sU,sL)=r,s
    if len(rU)>len(sU) or len(rL)>len(sL): return False
    m={}
    for R,S in ((rU,sU),(rL,sL)):
        for i in range(len(R)):
            rid=R[len(R)-1-i]; sid=S[len(S)-1-i]
            if rid==0: continue
            if rid in m and m[rid]!=sid: return False
            m[rid]=sid
    return True
def canon2(st):
    U,L=st;m={};nx=1
    for x in U+L:
        if x not in m: m[x]=nx;nx+=1
    return (tuple(m[x] for x in U),tuple(m[x] for x in L))
W=int(sys.argv[1]);k=int(sys.argv[2]);TOT=int(sys.argv[3]) if len(sys.argv)>3 else 0
start=(((),()),((),()));seen={start};stack=[start];bad=0;n=0
while stack:
    r,s=stack.pop()
    rn=rs(r,k);sn=es(s,W)
    for ev,t in sn.items():
        n+=1
        if ev not in rn or not covers(rn[ev],t): bad+=1; 
        else:
            p=(canon2(rn[ev]),canon2(t))
            # note: canonicalising r and t separately can break the joint id map; keep joint canon instead
            m={0:0};nx=1
            for x in rn[ev][0]+rn[ev][1]+t[0]+t[1]:
                if x not in m: m[x]=nx;nx+=1
            p=((tuple(m[x] for x in rn[ev][0]),tuple(m[x] for x in rn[ev][1])),(tuple(m[x] for x in t[0]),tuple(m[x] for x in t[1])))
            if p not in seen: seen.add(p);stack.append(p)
print("W",W,"k",k,"TOT",TOT,"pairs",len(seen),"transitions",n,"violations",bad)
cnt=0;seen={start};stack=[start]
while stack and cnt<4:
    r,s=stack.pop();rn=rs(r,k);sn=es(s,W)
    for ev,t in sn.items():
        if ev not in rn or not covers(rn[ev],t):
            cnt+=1;print("r",r,"s",s,"ev",ev,"-> r'",rn.get(ev),"t",t)
        else:
            m={};nx=1
            for x in rn[ev][0]+rn[ev][1]+t[0]+t[1]:
                if x not in m: m[x]=nx;nx+=1
            m[0]=0
            p=((tuple(m[x] for x in rn[ev][0]),tuple(m[x] for x in rn[ev][1])),(tuple(m[x] for x in t[0]),tuple(m[x] for x in t[1])))
            if p not in seen: seen.add(p);stack.append(p)
