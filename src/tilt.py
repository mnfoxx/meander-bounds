# tilted relaxation: edge weight theta^{dwidth}; rho(T_theta) is an upper bound per crossing for every theta>0
import sys, numpy as np, scipy.sparse as sp
from relax_py2 import canon, strip
from collections import Counter
def succ(state,k,TOT):
    U,L=state;out=[]
    for eu in (0,1):
        for el in (0,1):
            U2=list(U);L2=list(L);closed=[]
            if eu==0: closed.append(U2.pop() if U2 else 0)
            if el==0: closed.append(L2.pop() if L2 else 0)
            if len(closed)==2 and closed[0] and closed[0]==closed[1]: continue
            cur=0
            if len(closed)==2:
                a,b=closed
                if a and b: U2=[a if x==b else x for x in U2];L2=[a if x==b else x for x in L2];cur=a
            elif closed: cur=closed[0]
            if eu or el:
                if not cur: cur=max(U2+L2+[0])+1
                if eu: U2.append(cur)
                if el: L2.append(cur)
            U2,L2=strip(U2,L2)
            if TOT>0:
                while len(U2)+len(L2)>TOT:
                    if len(U2)>=len(L2): U2.pop(0)
                    else: L2.pop(0)
            else:
                while len(U2)>k: U2.pop(0)
                while len(L2)>k: L2.pop(0)
            U2,L2=strip(U2,L2);U2,L2=strip(U2,L2)
            out.append((canon(tuple(U2),tuple(L2)),(eu+el-1)*2))
    return out
