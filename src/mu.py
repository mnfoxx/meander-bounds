import sys, numpy as np
src=open('meanderphase.py').read().split('N=len(order)')[0]
for H in range(2,8):
    sys.argv=['x',str(H)]; ns={}; exec(src,ns)
    order=ns['order']; states=ns['states']; edges=ns['edges']; N=len(order)
    M=np.zeros((N,N))
    for s,ts in edges.items():
        for t in ts: M[states[s],states[t]]+=1
    mu=max(abs(np.linalg.eigvals(M)))
    print("D=%d states=%d mu_D=%.5f"%(H-1,N,mu))
