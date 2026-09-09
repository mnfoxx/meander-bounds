import sys
from itertools import product
def canon(R,L):
    m={};nxt=1
    def f(c):
        nonlocal nxt
        if c not in m: m[c]=nxt;nxt+=1
        return m[c]
    return (tuple(f(c) for c in R),tuple(f(c) for c in L))
def count(W,kmax):
    def step(st):
        R,L=st; outs=[]
        for er,el in product([1,0],repeat=2):
            if not er and not R: continue
            if not el and not L: continue
            if len(R)+len(L)+(1 if er else -1)+(1 if el else -1)>W: continue
            R2=list(R);L2=list(L);touched=[]
            if not er: touched.append(R2.pop())
            if not el: touched.append(L2.pop())
            if len(touched)==2 and touched[0]==touched[1]:
                if not R2 and not L2: outs.append('DONE')
                continue
            if not touched: cid=max(list(R2)+list(L2)+[0])+1
            else:
                cid=min(touched); R2=[cid if c in touched else c for c in R2]; L2=[cid if c in touched else c for c in L2]
            if er: R2.append(cid)
            if el: L2.append(cid)
            outs.append(canon(tuple(R2),tuple(L2)))
        return outs
    init=((),()); succ={}; order=[init]; i=0
    while i<len(order):
        s=order[i]; i+=1; succ[s]=step(s)
        for t in succ[s]:
            if t!='DONE' and t not in succ and t not in order[i:]: order.append(t)
    vec={init:1}; res=[]
    for m in range(1,2*kmax+1):
        nv={}; done=0
        for s,c in vec.items():
            for t in succ[s]:
                if t=='DONE': done+=c
                else: nv[t]=nv.get(t,0)+c
        vec=nv; res.append(done)
    return len(order),[res[2*k-1] for k in range(1,kmax+1)]
if __name__=="__main__":
  for W in (4,6,8):
      N,c=count(W,7); print("width<=%d states=%d:"%(W,N),c)
