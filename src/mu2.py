import sys, numpy as np
from scipy.sparse import lil_matrix, csr_matrix
src=open('meanderphase.py').read().split('N=len(order)')[0]
res={}
for H in range(2,10):
    sys.argv=['x',str(H)]; ns={}; exec(src,ns)
    order=ns['order']; states=ns['states']; edges=ns['edges']; N=len(order)
    M=lil_matrix((N,N))
    for s,ts in edges.items():
        for t in ts: M[states[s],states[t]]+=1
    M=csr_matrix(M); M2=M@M  # bipartite: iterate M^2
    v=np.ones(N); 
    for it in range(4000):
        w=M2@v; lam=np.linalg.norm(w); v=w/lam
    mu=np.sqrt((v@(M2@v))/(v@v)); res[H-1]=mu
    print("D=%d states=%d mu_D=%.6f"%(H-1,N,mu),flush=True)
