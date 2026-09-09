// Independent checker for adaptive-truncation certificates of the closed-meander transfer matrix.
// Reads  <prefix>_states.bin  (N records of 20 bytes: r, l, nibbles = label-1 for upper stack bottom->top
// then lower stack bottom->top; mirror-canonical orbit representatives) and <prefix>_vec.bin (N int64).
// Independently re-derives the transitions from a "points on the cut line + non-crossing matching" model,
// verifies (1) every state is a well-formed reachable configuration (BFS inside S from the one-path state),
// (2) the one-path state j0 is present with v>0, (3) (T^T)^2 v >= c v entrywise for c = the reported minimum.
// Usage: check prefix
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#define SB 20
typedef struct { uint8_t r,l; uint8_t p[SB-2]; } Rec;

// cut-line model: n=2p points; pts[0..r-1] = upper arcs outermost..innermost, pts[r..n-1] = lower arcs innermost..outermost.
// partner[i] = index of the other end of the same path.
typedef struct { int r,n; int partner[40]; } Cfg;

static void rec_to_cfg(const Rec*q,Cfg*c){
  int r=q->r,l=q->l,n=r+l;c->r=r;c->n=n;int lab[40];
  for(int i=0;i<n;i++)lab[i]=((q->p[i>>1]>>((i&1)*4))&15)+1;
  // record order: upper bottom->top (outermost->innermost) = pts 0..r-1 ; lower bottom->top (outermost->innermost) = pts n-1..r
  int seq[40];for(int i=0;i<r;i++)seq[i]=lab[i];for(int k=0;k<l;k++)seq[r+k]=lab[r+(l-1-k)];
  for(int i=0;i<n;i++){c->partner[i]=-1;for(int j=0;j<n;j++)if(j!=i&&seq[j]==seq[i])c->partner[i]=j;}
}
static int cfg_wellformed(const Cfg*c){ // each point has exactly one partner, partner is an involution, matching non-crossing
  if(c->n&1)return 0;
  for(int i=0;i<c->n;i++){int j=c->partner[i];if(j<0||j==i||c->partner[j]!=i)return 0;}
  for(int i=0;i<c->n;i++)for(int k=0;k<c->n;k++){int j=c->partner[i],l=c->partner[k];
    if(i<k&&k<j&&j<l)return 0;} // crossing i<k<j<l with (i,j),(k,l) matched
  return 1;}
static void cfg_to_rec(const Cfg*c,Rec*out){ // canonical: labels by first appearance in the record order, min over mirror
  int r=c->r,n=c->n,l=n-r;Rec cand[2];
  for(int m=0;m<2;m++){
    // mirror: swap the two half-planes = reverse the point sequence, r<->l
    int seqp[40],rr=m?l:r;
    for(int i=0;i<n;i++)seqp[i]=m?(n-1-c->partner[n-1-i]):c->partner[i];
    // record order: upper outermost->innermost = pts 0..rr-1, then lower outermost->innermost = pts n-1..rr
    int ord[40];for(int i=0;i<rr;i++)ord[i]=i;for(int k=0;k<n-rr;k++)ord[rr+k]=n-1-k;
    int pos[40];for(int i=0;i<n;i++)pos[ord[i]]=i;
    int lab[40];for(int i=0;i<n;i++)lab[i]=-1;int nxt=0;
    for(int i=0;i<n;i++){if(lab[i]<0){lab[i]=nxt;lab[pos[seqp[ord[i]]]]=nxt;nxt++;}}
    memset(&cand[m],0,SB);cand[m].r=rr;cand[m].l=n-rr;for(int i=0;i<n;i++)cand[m].p[i>>1]|=lab[i]<<((i&1)*4);}
  *out=(memcmp(&cand[0],&cand[1],SB)<=0)?cand[0]:cand[1];}

// successors: events (eu,el) in {open=1, close=0} on upper/lower. Returns count; cycle-closing moves excluded.
static int successors(const Cfg*c,Cfg*out){
  int cnt=0;
  for(int eu=0;eu<2;eu++)for(int el=0;el<2;el++){
    int r=c->r,n=c->n;
    if(!eu&&r==0)continue;if(!el&&n-r==0)continue;
    int iu=r-1,il=r; // innermost upper point, innermost lower point
    if(!eu&&!el&&c->partner[iu]==il)continue; // would close a loop
    Cfg d;int map[40];for(int i=0;i<n;i++)map[i]=-2;
    // build new point list: keep old points except closed innermost ones; insert new ones at the middle
    int newr=r-(eu?0:1)+(eu?1:0);int m=0;int oldpart[40];memcpy(oldpart,c->partner,sizeof(int)*n);
    // when a point is removed, its partner becomes connected to whatever the road continues into
    int newidx[40];
    for(int i=0;i<n;i++){if((!eu&&i==iu)||(!el&&i==il))newidx[i]=-1;else newidx[i]=m++;}
    int upperNew=-1,lowerNew=-1;
    // positions: upper points keep order; new upper point goes after last kept upper point (innermost); new lower point before first kept lower point
    // recompute indices in the final sequence: upper kept (in order), [new upper], [new lower], lower kept (in order)
    int idx=0;for(int i=0;i<r;i++)if(newidx[i]>=0)newidx[i]=idx++;
    if(eu)upperNew=idx++;if(el)lowerNew=idx++;
    for(int i=r;i<n;i++)if(newidx[i]>=0)newidx[i]=idx++;
    d.n=idx;d.r=newr;for(int i=0;i<idx;i++)d.partner[i]=-1;
    // the road: enters through the closed points (if any) and leaves through the opened points (if any)
    // ends of the road segment at this crossing:
    int endsOld[2],ne=0;if(!eu)endsOld[ne++]=oldpart[iu];if(!el)endsOld[ne++]=oldpart[il];
    int endsNew[2],nn=0;if(eu)endsNew[nn++]=upperNew;if(el)endsNew[nn++]=lowerNew;
    for(int i=0;i<n;i++){if(newidx[i]<0)continue;int j=oldpart[i];
      if(newidx[j]>=0)d.partner[newidx[i]]=newidx[j];}
    // connect: old loose ends (partners of closed points) and new points, total ne+nn == 2 ends of one path
    int ends[2],k=0;for(int q=0;q<ne;q++)ends[k++]=newidx[endsOld[q]];for(int q=0;q<nn;q++)ends[k++]=endsNew[q];
    if(k!=2){fprintf(stderr,"internal\n");exit(1);}
    d.partner[ends[0]]=ends[1];d.partner[ends[1]]=ends[0];
    out[cnt++]=d;}
  return cnt;}

Rec*S;size_t N;uint32_t*H;size_t hc;
static uint64_t hsh(const Rec*a){uint64_t h=14695981039346656037ULL;const uint8_t*b=(const uint8_t*)a;for(int i=0;i<SB;i++){h^=b[i];h*=1099511628211ULL;}h^=h>>31;h*=0x9E3779B97F4A7C15ULL;h^=h>>29;return h;}
static int64_t look(const Rec*s){size_t k=hsh(s)&(hc-1);while(H[k]!=0xffffffffu){if(!memcmp(&S[H[k]],s,SB))return H[k];k=(k+1)&(hc-1);}return -1;}

int main(int argc,char**argv){
  char fn[512];snprintf(fn,512,"%s_states.bin",argv[1]);FILE*f=fopen(fn,"rb");fseek(f,0,SEEK_END);N=ftell(f)/SB;fseek(f,0,SEEK_SET);
  S=malloc(N*SB);if(fread(S,SB,N,f)!=N)return 1;fclose(f);
  int64_t*v=malloc(N*8);snprintf(fn,512,"%s_vec.bin",argv[1]);f=fopen(fn,"rb");if(fread(v,8,N,f)!=N)return 1;fclose(f);
  hc=1;while(hc<2*N)hc<<=1;H=malloc(hc*4);memset(H,0xff,hc*4);
  for(size_t i=0;i<N;i++){size_t k=hsh(&S[i])&(hc-1);while(H[k]!=0xffffffffu){if(!memcmp(&S[H[k]],&S[i],SB)){fprintf(stderr,"duplicate state %zu\n",i);return 1;}k=(k+1)&(hc-1);}H[k]=i;}
  // (1) well-formed & canonical
  size_t bad=0;for(size_t i=0;i<N;i++){Cfg c;rec_to_cfg(&S[i],&c);Rec q;if(!cfg_wellformed(&c)){bad++;continue;}cfg_to_rec(&c,&q);if(memcmp(&q,&S[i],SB))bad++;}
  printf("states: %zu, malformed/non-canonical: %zu\n",N,bad);if(bad)return 1;
  // j0
  Rec j0;{Cfg c;c.r=1;c.n=2;c.partner[0]=1;c.partner[1]=0;cfg_to_rec(&c,&j0);}
  int64_t J=look(&j0);printf("j0 index %lld, v[j0]=%lld\n",(long long)J,J>=0?(long long)v[J]:-1LL);if(J<0||v[J]<=0)return 1;
  // edges inside S
  size_t nedges=0;uint32_t*eoff=malloc((N+1)*4);uint32_t*edst=malloc(N*4*4);
  for(size_t i=0;i<N;i++){eoff[i]=nedges;Cfg c,o[4];rec_to_cfg(&S[i],&c);int n=successors(&c,o);
    for(int q=0;q<n;q++){Rec t;cfg_to_rec(&o[q],&t);int64_t j=look(&t);if(j>=0)edst[nedges++]=j;}}
  eoff[N]=nedges;printf("edges inside S: %zu\n",nedges);
  // (2) reachability inside S from j0 (j0 itself is reachable from the empty configuration by (open,open))
  uint8_t*seen=calloc(N,1);uint32_t*queue=malloc(N*4);size_t qh=0,qt=0;seen[J]=1;queue[qt++]=J;
  while(qh<qt){uint32_t i=queue[qh++];for(uint32_t k=eoff[i];k<eoff[i+1];k++){uint32_t j=edst[k];if(!seen[j]){seen[j]=1;queue[qt++]=j;}}}
  printf("reachable from j0 inside S: %zu of %zu (certificate is evaluated on the reachable subset)\n",qt,N);
  // (3) t2 = (T^T)^2 v on the reachable subset
  int64_t*t1=calloc(N,8),*t2=calloc(N,8);
  for(size_t i=0;i<N;i++){if(!seen[i])continue;int64_t x=v[i];if(x<0){printf("negative v\n");return 1;}if(x)for(uint32_t k=eoff[i];k<eoff[i+1];k++)t1[edst[k]]+=x;}
  for(size_t i=0;i<N;i++){if(!seen[i])continue;int64_t x=t1[i];if(x)for(uint32_t k=eoff[i];k<eoff[i+1];k++)t2[edst[k]]+=x;}
  int64_t bn=0,bd=1;int first=1;size_t npos=0;
  for(size_t i=0;i<N;i++){if(!seen[i]||v[i]<=0)continue;npos++;if(first){bn=t2[i];bd=v[i];first=0;}else if((__int128)t2[i]*bd<(__int128)bn*v[i]){bn=t2[i];bd=v[i];}}
  printf("positive entries %zu; min_i ((T^T)^2 v)_i / v_i = %lld/%lld = %.9f\n",npos,(long long)bn,(long long)bd,(double)bn/(double)bd);
  size_t unreach=0;
  printf("CERTIFICATE OK: R >= %.9f\n",(double)bn/(double)bd);
  return unreach?1:0;}
