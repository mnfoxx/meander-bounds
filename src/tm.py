# Transfer matrix for depth<=1 semi-meanders, scanning pile top->bottom.
# State: (Rstack, Lstack, tail_used) where stacks are tuples of component ids (bottom..top),
# component ids canonicalised; id 0 = root component (contains stamp 1, one loose end).
# tail component (contains stamp n, one loose end) gets id 'T'. Others have 2 loose ends.
from itertools import product
from sympy import Matrix, symbols, simplify, eye, factor, cancel, series, Rational, fibonacci, sqrt

MAXH=2
def canon(R,L,tail):
    # relabel non-root, non-tail components in order of first appearance
    m={}; nxt=1
    def f(c):
        nonlocal nxt
        if c in (0,'T'): return c
        if c not in m: m[c]=nxt; nxt+=1
        return m[c]
    return (tuple(f(c) for c in R), tuple(f(c) for c in L), tail)

def step(state):
    """yield successor states after processing one position (one stamp)."""
    R,L,tail=state
    outs=[]
    # each side: 'open', 'close', 'absent'
    for er,el in product(['open','close','absent'],repeat=2):
        if er=='absent' and el=='absent': continue
        if er=='absent' and tail: continue   # tail slot may be used once
        if el=='absent' and tail: continue
        if er=='close' and not R: continue
        if el=='close' and not L: continue
        if er=='open' and len(R)>=MAXH: continue
        if el=='open' and len(L)>=MAXH: continue
        R2=list(R); L2=list(L); ntail=tail
        # components touched by this stamp
        touched=[]
        if er=='close': touched.append(R2.pop())
        if el=='close': touched.append(L2.pop())
        if er=='absent' or el=='absent': ntail=True
        # merge rule: the stamp joins the touched components (and tail if absent)
        if len(touched)==2:
            a,b=touched
            if a==b: continue  # cycle
        # new id for the stamp's component
        newid=max([c for c in R2+L2 if c not in (0,'T')]+[0])+1
        # determine resulting component id
        ids=set(touched)
        if er=='absent' or el=='absent': ids.add('T')
        if 0 in ids and 'T' in ids:
            # root and tail merge: path complete; only valid if nothing else open and no opens now
            if R2 or L2 or er=='open' or el=='open': continue
            outs.append('DONE'); continue
        if 0 in ids: cid=0
        elif 'T' in ids: cid='T'
        elif ids: cid=min(ids)
        else: cid=newid
        # relabel touched ids to cid
        R2=[cid if c in ids else c for c in R2]
        L2=[cid if c in ids else c for c in L2]
        # opens push cid
        if er=='open': R2.append(cid)
        if el=='open': L2.append(cid)
        # a component with a single loose end and no root/tail role is impossible unless it is root/tail:
        # if cid is fresh (no touched, no absent) and no opens -> isolated stamp: invalid
        if not ids and er!='open' and el!='open': continue
        # a non-root/tail component must retain exactly 2 loose ends; root/tail exactly 1
        cnt={}
        for c in R2+L2: cnt[c]=cnt.get(c,0)+1
        ok=True
        for c,k in cnt.items():
            need=1 if c in (0,'T') else 2
            if k!=need: ok=False
        if cnt.get(0,0)!=1: ok=False
        if ntail and cnt.get('T',0)!=1: ok=False
        if not ok: continue
        outs.append(canon(tuple(R2),tuple(L2),ntail))
    return outs

# initial: position 1 processed: R open with root, L absent
init=((0,),(),False)
states={init:0}; order=[init]; edges={}
i=0
while i<len(order):
    s=order[i]; i+=1
    for t in step(s):
        if t=='DONE': edges.setdefault(s,[]).append('DONE'); continue
        if t not in states: states[t]=len(order); order.append(t)
        edges.setdefault(s,[]).append(t)
N=len(order); print("states:",N)
for s in order: print(states[s],s,'->',[ 'DONE' if t=='DONE' else states[t] for t in edges.get(s,[])])
T=Matrix(N,N,lambda i,j:0); d=Matrix(N,1,lambda i,j:0)
for s,ts in edges.items():
    for t in ts:
        if t=='DONE': d[states[s]]+=1
        else: T[states[s],states[t]]+=1
# a(n) = e0^T T^(n-2) d   (n-1 transitions: n-2 internal steps then final DONE step)
v=Matrix(1,N,lambda i,j:1 if j==0 else 0)
cnt=[]
w=v
for n in range(2,30):
    cnt.append((w*d)[0]); w=w*T
print("counts n=2..:",cnt)
x=symbols('x')
G=(v*(eye(N)-x*T).inv()*d)[0]*x**2
G=cancel(G); print("GF:",factor(G))
n=symbols('n')
def formula(n): return 2*fibonacci(n+2)-(3*2**(n//2) if n%2==0 else 2**((n+3)//2))
print("formula:",[formula(k) for k in range(3,30)])
print("match:",[formula(k) for k in range(3,30)]==cnt[1:])
