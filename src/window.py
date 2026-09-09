# Window certificates for lower bounds on the meander constant (prototype).
# Pattern = (U, L, fU, fL): U, L = top-k window entries (outer->inner), 0 = '?' (partner outside window),
# labels 1.. appear exactly twice.  fU = 1 iff there are deep (forgotten) arcs below the U window.
# v(exact state) = g(pattern) * theta^(#deep arcs).  We compute the best c with  sum_{certified predecessors} v(s) >= c v(t).
import sys, itertools
from collections import defaultdict

def canon(U,L,fU,fL,mirror=True):
    def relab(U,L):
        m={0:0};nx=1;out=[]
        for x in U+L:
            if x not in m: m[x]=nx;nx+=1
        return (tuple(m[x] for x in U),tuple(m[x] for x in L))
    a=relab(U,L)+(fU,fL)
    if not mirror: return a
    b=relab(L,U)+(fL,fU)
    return min(a,b)

def valid(U,L):
    seq=list(U)+list(reversed(L))   # cut-line order
    pos=defaultdict(list)
    for i,x in enumerate(seq):
        if x: pos[x].append(i)
    ch=[]
    for x,p in pos.items():
        if len(p)!=2: return False
        ch.append(tuple(p))
    for (a,b) in ch:
        for (c,d) in ch:
            if a<c<b<d: return False
    for i,x in enumerate(seq):
        if x==0:
            for (a,b) in ch:
                if a<i<b: return False
    return True

def first_q_from_outer(S):   # index of first '?' scanning outer->inner (index 0 is outer)
    for i,x in enumerate(S):
        if x==0: return i
    return -1
def first_q_from_inner(S):
    for i in range(len(S)-1,-1,-1):
        if S[i]==0: return i
    return -1

def refill_options(U,L,side,k):
    """side 0: a deep entry enters at the outer end of U (U currently has k-1 entries). Returns list of (U',L') options."""
    U=list(U);L=list(L);opts=[]
    new=max(U+L+[0])+1
    if side==0:
        # (i) '?'
        opts.append(([0]+U,L))
        # (ii) partner = first '?' in U from outer end
        i=first_q_from_outer(U)
        if i>=0:
            U2=U[:];U2[i]=new;opts.append(([new]+U2,L))
        # (iii) partner = first '?' of L from the inner end, requires no '?' in U
        if i<0:
            j=first_q_from_inner(L)
            if j>=0:
                L2=L[:];L2[j]=new;opts.append(([new]+U,L2))
    else:
        opts.append((U,[0]+L))
        i=first_q_from_outer(L)
        if i>=0:
            L2=L[:];L2[i]=new;opts.append((U,[new]+L2))
        if i<0:
            j=first_q_from_inner(U)
            if j>=0:
                U2=U[:];U2[j]=new;opts.append((U2,[new]+L))
    return opts

def push_out(S):
    """window overflow: drop the outer entry; its partner (if labelled) becomes '?'. returns new list"""
    x=S[0];S=S[1:]
    return S,x

def preds(t,k):
    """certified predecessors of pattern t: list of (list of option patterns, delta_d). The certificate uses min over options."""
    U,L,fU,fL=t;U=list(U);L=list(L);res=[]
    # ---- (A) reverse (open,open)
    if U and L and U[-1]!=0 and U[-1]==L[-1]:
        U1=U[:-1];L1=L[:-1];dd=-(fU+fL)
        opts=[(U1,L1)]
        if fU: opts=[o for (u,l) in opts for o in refill_options(u,l,0,k)]
        if fL:
            opts2=[]
            for (u,l) in opts: opts2+=refill_options(u,l,1,k)
            opts=opts2
            if fU and 0 not in U1 and 0 not in L1:
                new=max(U1+L1+[0])+1
                opts.append(([new]+U1,[new]+L1))
        res.append((opts,dd,'?' if fU else '0','?' if fL else '0'))
    # ---- (B) reverse (close U, open L): t's L-top moves to U-top
    if L:
        c=L[-1];L1=L[:-1];dd=0
        Ufull=len(U)==k
        U1=U+[c]
        # the new U top: label c (if c!=0 it stays paired with c's other end); if c==0 it's '?'
        opts=[(U1,L1)]
        # option: c==0 and its partner is exactly the entering deep L entry (handled below via refill: pair new U top with entering)
        if Ufull:
            # push out outer U entry
            opts2=[]
            for (u,l) in opts:
                x=u[0];u2=u[1:]
                if x:
                    u2=[0 if y==x else y for y in u2];l=[0 if y==x else y for y in l]
                opts2.append((u2,l))
            opts=opts2;dd+=1
        if fL:
            dd-=1
            opts2=[]
            for (u,l) in opts:
                opts2+=refill_options(u,l,1,k)
                if c==0:
                    # new U top (u[-1]==0) paired with the entering L entry
                    if u and u[-1]==0:
                        new=max(u+l+[0])+1
                        u2=u[:];u2[-1]=new
                        # entering entry at outer end of L; chord (U-top, L-outer) encloses everything: requires all internal
                        if 0 not in u2 and 0 not in l:
                            opts2.append((u2,[new]+l))
            opts=opts2
        res.append((opts,dd,'1' if Ufull else '0','?' if fL else '0'))
    # ---- (C) reverse (open U, close L): symmetric
    if U:
        c=U[-1];U1=U[:-1];dd=0
        Lfull=len(L)==k
        L1=L+[c]
        opts=[(U1,L1)]
        if Lfull:
            opts2=[]
            for (u,l) in opts:
                x=l[0];l2=l[1:]
                if x:
                    l2=[0 if y==x else y for y in l2];u=[0 if y==x else y for y in u]
                opts2.append((u,l2))
            opts=opts2;dd+=1
        if fU:
            dd-=1
            opts2=[]
            for (u,l) in opts:
                opts2+=refill_options(u,l,0,k)
                if c==0 and l and l[-1]==0:
                    new=max(u+l+[0])+1
                    l2=l[:];l2[-1]=new
                    if 0 not in u and 0 not in l2:
                        opts2.append(([new]+u,l2))
            opts=opts2
        res.append((opts,dd,'?' if fU else '0','1' if Lfull else '0'))
    # ---- (D) reverse (close,close): split a visible pair (x,y) into a~x, b~y or a~y, b~x
    labels=set(x for x in U+L if x)
    Ufull=len(U)==k;Lfull=len(L)==k
    for lab in labels:
        for orient in (0,1):
            new=max(U+L)+1
            # positions of lab: replace one occurrence by 'new' (partner of a), the other keeps lab (partner of b)
            occ=[(0,i) for i,x in enumerate(U) if x==lab]+[(1,i) for i,x in enumerate(L) if x==lab]
            (s1,i1),(s2,i2)=occ
            if orient: (s1,i1),(s2,i2)=(s2,i2),(s1,i1)
            U2=U[:];L2=L[:]
            if s1==0: U2[i1]=new
            else: L2[i1]=new
            # a on U-top with label new ; b on L-top with label lab
            U2=U2+[new];L2=L2+[lab]
            dd=0;fU2=fU;fL2=fL
            if Ufull:
                x=U2[0];U2=U2[1:]
                if x: U2=[0 if y==x else y for y in U2];L2=[0 if y==x else y for y in L2]
                dd+=1;fU2=1
            if Lfull:
                x=L2[0];L2=L2[1:]
                if x: L2=[0 if y==x else y for y in L2];U2=[0 if y==x else y for y in U2]
                dd+=1;fL2=1
            if valid(U2,L2):
                res.append(([(U2,L2)],dd,'1' if Ufull else '0','1' if Lfull else '0'))
    return res

def enumerate_patterns(k):
    """all valid patterns: sequences over atoms '?' / chord-groups, split into U|L with each side <= k, flags consistent."""
    # generate complete non-crossing matchings on n points as chord lists
    from functools import lru_cache
    @lru_cache(None)
    def matchings(n):
        if n==0: return [()]
        res=[]
        for i in range(1,n,2):
            for m1 in matchings(i-1):
                for m2 in matchings(n-i-1):
                    res.append(((0,i),)+tuple((a+1,b+1) for a,b in m1)+tuple((a+i+1,b+i+1) for a,b in m2))
        return res
    seqs={0:[[]]}   # length -> list of sequences (as list of entries with labels)
    def gen(n):
        if n in seqs: return seqs[n]
        out=[]
        for s in gen(n-1): out.append([0]+s)
        for j in range(1,n//2+1):
            for m in matchings(2*j-2):
                # chord group of length 2j: outer chord (0,2j-1) enclosing matching m on inner points 1..2j-2
                grp=[None]*(2*j)
                grp[0]=grp[2*j-1]='G'
                for (a,b) in m: grp[a+1]=grp[b+1]=(a,b)
                for s in gen(n-2*j): out.append(grp+s)
        seqs[n]=out;return out
    pats=set()
    for n in range(0,2*k+1):
        for s in gen(n):
            # relabel: distinct objects -> labels
            nx=1;lbl=[0]*n;i=0
            while i<n:
                if s[i]==0: lbl[i]=0;i+=1;continue
                # find group extent: s[i]=='G' begins a group; its end is the matching 'G' after the enclosed matching
                # groups were placed contiguously: scan forward to find the closing 'G' at the same nesting
                j=i+1
                while s[j]!='G': j+=1
                # chord (i,j) plus inner matching entries encoded as (a,b) relative to i+1
                lbl[i]=lbl[j]=nx;nx+=1
                inner=s[i+1:j]
                seen={}
                for q,x in enumerate(inner):
                    if x not in seen: seen[x]=nx;nx+=1
                    lbl[i+1+q]=seen[x]
                i=j+1
            for r in range(0,n+1):
                if r>k or n-r>k: continue
                U=tuple(lbl[:r]);L=tuple(reversed(lbl[r:]))
                for fU in (0,1):
                    for fL in (0,1):
                        if fU and r<k: continue
                        if fL and n-r<k: continue
                        if (0 in U or 0 in L) and not (fU or fL): continue
                        pats.add(canon(U,L,fU,fL))
    return sorted(pats)

if __name__=="__main__":
    k=int(sys.argv[1]);thetas=[float(x) for x in sys.argv[2].split(',')];iters=int(sys.argv[3]) if len(sys.argv)>3 else 300
    pats=enumerate_patterns(k);idx={p:i for i,p in enumerate(pats)};n=len(pats)
    print("k",k,"patterns",n,flush=True)
    # precompute predecessor structure: for each t: list of (list of option indices, dd)
    import math
    struct=[]
    missing=0
    for t in pats:
        entry=[]
        for opts,dd,sU,sL in preds(t,k):
            ids=[]
            for (u,l) in opts:
                cands=[(a,b) for a in ([0,1] if sU=='?' else [int(sU)]) for b in ([0,1] if sL=='?' else [int(sL)])]
                for a,b in cands:
                    if not valid(u,l): continue
                    if (0 in u or 0 in l) and not (a or b): continue
                    if a and len(u)<k: continue
                    if b and len(l)<k: continue
                    key=canon(tuple(u),tuple(l),a,b)
                    if key in idx: ids.append(idx[key])
                    else: missing+=1
            if ids: entry.append((ids,dd))
        struct.append(entry)
    print("missing option patterns:",missing,flush=True)
    j0=idx[canon((1,),(1,),0,0)]
    flagged=[1 if (p[2] or p[3]) else 0 for p in pats]
    def Phi(g,theta):
        ng=[0.0]*n
        for i in range(n):
            if theta==0 and flagged[i]: continue
            s=0.0
            for ids,dd in struct[i]:
                if theta==0:
                    if dd==0: s+=min(g[j] for j in ids)
                else: s+=min(g[j] for j in ids)*theta**dd
            ng[i]=s
        return ng
    def feasible(c,theta,iters):
        g=[0.0 if (theta==0 and flagged[i]) else 1.0 for i in range(n)]
        for it in range(iters):
            pg=Phi(g,theta)
            g2=[min(g[i],pg[i]/c) for i in range(n)]
            if g2[j0]<1e-9: return None
            g=g2
        pg=Phi(g,theta)
        ratio=min(pg[i]/g[i] for i in range(n) if g[i]>1e-12)
        return ratio,g
    for theta in thetas:
        # bisection on c with the decreasing (greatest-fixed-point) iteration; certify with the final g
        lo,hi=1.0,3.6;best=None
        for _ in range(14):
            c=(lo+hi)/2
            r=feasible(c,theta,iters)
            if r is None: hi=c
            else:
                lo=c
                if best is None or r[0]>best[0]: best=r
        ratio,g=best
        npos=sum(1 for x in g if x>1e-12)
        print(f"k={k} theta={theta} certified c={ratio:.5f} c^2={ratio*ratio:.4f} (feasible up to {lo:.4f}) positive={npos}/{n}",flush=True)
