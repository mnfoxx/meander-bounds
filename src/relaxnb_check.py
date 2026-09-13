# Mechanical check of Proposition 8.5 (paper v2): every word accepted by the restricted relaxation R''_m
# (no close on an empty partial stack, no (close,close) popping a '?') is a valid cycle-free walk of the
# exact automaton. Enumerates all accepted words up to length N.   Usage: python3 relaxnb_check.py m N
import sys
from collections import Counter
sys.setrecursionlimit(10000)
def strip(U,L):
    c=Counter([x for x in U+L if x])
    U=[x if (x and c[x]==2) else 0 for x in U];L=[x if (x and c[x]==2) else 0 for x in L]
    while U and U[0]==0: U.pop(0)
    while L and L[0]==0: L.pop(0)
    return U,L
def rstep(state,eu,el,TOT):
    U,L=state;U2=list(U);L2=list(L);closed=[]
    if (eu==0 and not U2) or (el==0 and not L2): return None
    if eu==0: closed.append(U2.pop())
    if el==0: closed.append(L2.pop())
    if len(closed)==2 and (closed[0]==0 or closed[1]==0): return None
    if len(closed)==2 and closed[0]==closed[1]: return None
    cur=0
    if len(closed)==2:
        a,b=closed;U2=[a if x==b else x for x in U2];L2=[a if x==b else x for x in L2];cur=a
    elif closed: cur=closed[0]
    if eu or el:
        if not cur: cur=max(U2+L2+[0])+1
        if eu: U2.append(cur)
        if el: L2.append(cur)
    U2,L2=strip(U2,L2)
    while len(U2)+len(L2)>TOT:
        if len(U2)>=len(L2): U2.pop(0)
        else: L2.pop(0)
    U2,L2=strip(U2,L2);U2,L2=strip(U2,L2)
    return (tuple(U2),tuple(L2))
def estep(state,eu,el):
    U,L=state;U2=list(U);L2=list(L);closed=[]
    if eu==0:
        if not U2: return None
        closed.append(U2.pop())
    if el==0:
        if not L2: return None
        closed.append(L2.pop())
    if len(closed)==2 and closed[0]==closed[1]: return None
    if len(closed)==2:
        a,b=closed;U2=[a if x==b else x for x in U2];L2=[a if x==b else x for x in L2];cur=a
    elif closed: cur=closed[0]
    else: cur=max(U2+L2+[0])+1
    if eu: U2.append(cur)
    if el: L2.append(cur)
    return (tuple(U2),tuple(L2))
m=int(sys.argv[1]);N=int(sys.argv[2]);bad=0;count=[0]*(N+1);A=[0]*(N+1)
def dfs(rs,es,n):
    global bad
    count[n]+=1
    if n==N: return
    for eu in (0,1):
        for el in (0,1):
            r2=rstep(rs,eu,el,m)
            if r2 is None: continue
            e2=estep(es,eu,el)
            if e2 is None: bad+=1;continue
            dfs(r2,e2,n+1)
def dfs2(es,n):
    A[n]+=1
    if n==N: return
    for eu in (0,1):
        for el in (0,1):
            e2=estep(es,eu,el)
            if e2 is not None: dfs2(e2,n+1)
dfs(((),()),((),()),0);dfs2(((),()),0)
print(f"m={m} N={N}: restricted words that are not valid cycle-free exact walks: {bad}")
print("restricted counts by length:",count)
print("cycle-free counts A_n       :",A)
