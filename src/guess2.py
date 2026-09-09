from sympy import Matrix, Rational
from guess import poly_fit
def lin_rec(seq,maxorder=7):
    N=len(seq)
    for r in range(1,maxorder+1):
        if N-r < r+3: break
        rows=[[Rational(seq[n-i]) for i in range(r+1)] for n in range(r,N)]
        ns=Matrix(rows).nullspace()
        if ns:
            v=ns[0]; v=v/v[0]; return r,[x for x in v]
    return None
hair={ # hairpins row per n, n=3..18 index0=0 hairpins
3:[0,1,1],4:[0,1,2,1],5:[0,1,5,3,1],6:[0,1,6,12,4,1],7:[0,1,11,26,22,5,1],8:[0,1,14,49,68,35,6,1],
9:[0,1,19,99,183,143,51,7,1],10:[0,1,22,158,402,486,258,70,8,1],11:[0,1,31,256,896,1404,1096,424,92,9,1],
12:[0,1,34,371,1646,3456,3774,2140,648,117,10,1],13:[0,1,45,561,3099,8246,11588,8914,3826,941,145,11,1],
14:[0,1,50,750,5072,17045,30628,31470,18428,6335,1310,176,12,1],
15:[0,1,57,1100,8626,34989,78057,100282,76660,35139,9945,1766,210,13,1],
16:[0,1,64,1401,13066,63931,174738,280955,274860,165391,61998,14891,2316,247,14,1],
17:[0,1,79,1901,20791,118773,386135,751281,904067,685259,331357,103749,21521,2971,287,15,1],
18:[0,1,84,2340,29808,198732,767446,1796420,2640786,2495502,1528826,615210,165102,30124,3738,330,16,1]}
for k in range(2,8):
    s=[hair[n][n-k] for n in range(max(3,k+2),19)]
    print(f"hairpins=n-{k}:",s)
    print("   poly deg all/even/odd:",poly_fit(s),poly_fit(s[0::2]),poly_fit(s[1::2]))
depth1=[2,2,8,16,34,60,112,190,336,560,962,1588,2680,4398,7336,11992]
print("depth1 linrec:",lin_rec(depth1))
depth12=[a+b for a,b in zip(depth1,[0,0,0,6,30,94,289,732,1935,4547,11226,25468,60557,134802,313133,689300])]
print("depth<=2 linrec:",lin_rec(depth12))
print("ratios depth1:",[round(depth1[i+1]/depth1[i],4) for i in range(8,15)])
