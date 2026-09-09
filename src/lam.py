import sys, numpy as np
src=open('tmD.py').read().split('import numpy')[0].replace('print("states:",N)','')
for H in [6]:
    sys.argv=['x',str(H)]; ns={}; exec(src,ns)
    order=ns['order']; edges=ns['edges']
    pre=[s for s in order if not s[2]]; idx={s:i for i,s in enumerate(pre)}
    M=np.zeros((len(pre),len(pre)))
    for s in pre:
        for t in edges.get(s,[]):
            if t!='DONE' and t in idx: M[idx[s],idx[t]]+=1
    print("D=%d pre-q states=%d lambda_D=%.5f"%(H-1,len(pre),max(abs(np.linalg.eigvals(M)))))
