// Upper bound on the meander constant R from the label-forgetting relaxation (per-side memory k).
// State: upper/lower stacks of labelled arcs (at most k per side; deeper arcs forgotten). Every closed meander
// is a walk of this automaton, so rho(T)^2 >= R.  Certificate: integer u>0 with T^2 u <= c u entrywise => rho(T^2) <= c.
// relaxnb_big: the RESTRICTED relaxation R''_m (no close on an empty partial stack, no (close,close) popping a '?').
// Its words are cycle-free exact walks, so its Perron root is a certified LOWER bound on the prefix constant gamma.
// Usage: relaxnb_big k iters total LOG2HASH STATECAP   e.g.  relaxnb_big 40 600 28 29 320000000
//   LOG2HASH: hash table has 2^LOG2HASH slots (need > 1.5 x states); STATECAP: preallocated state capacity.
// Large-machine version: 18-byte states (up to 32 entries), 64-bit edge offsets, two-pass edge build (no realloc
// spikes), OpenMP power iteration. Compile: gcc -O2 -fopenmp -o relax_big relax_big.c -lm
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#define MAXW 32
#define SB 18
typedef struct { uint8_t r,l; uint8_t p[16]; } St; // 18 bytes, nibble-packed labels 0..15 (0 = forgotten); up to 32 entries
int K,TOT;
static uint64_t hsh(const St*a){uint64_t h=1469598103934665603ULL;const uint8_t*b=(const uint8_t*)a;for(int i=0;i<SB;i++){h^=b[i];h*=1099511628211ULL;}h^=h>>29;h*=0x9E3779B97F4A7C15ULL;h^=h>>32;return h;}
static int eq(const St*a,const St*b){return memcmp(a,b,SB)==0;}
static inline int getl(const St*s,int i){return (s->p[i>>1]>>((i&1)*4))&15;}
static void pack(St*s,const uint8_t*lab,int n){memset(s->p,0,16);for(int i=0;i<n;i++)s->p[i>>1]|=lab[i]<<((i&1)*4);}
static void canon_lab(uint8_t*lab,int n){uint8_t map[64];memset(map,0,64);int nxt=1;for(int i=0;i<n;i++){uint8_t c=lab[i];if(!c)continue;if(!map[c])map[c]=nxt++;lab[i]=map[c];}}
static void make_state(St*t,const uint8_t*R,int nr,const uint8_t*L,int nl){
  uint8_t a[MAXW+2],b[MAXW+2];memcpy(a,R,nr);memcpy(a+nr,L,nl);canon_lab(a,nr+nl);memcpy(b,L,nl);memcpy(b+nl,R,nr);canon_lab(b,nr+nl);
  St t1,t2;memset(&t1,0,SB);memset(&t2,0,SB);t1.r=nr;t1.l=nl;pack(&t1,a,nr+nl);t2.r=nl;t2.l=nr;pack(&t2,b,nr+nl);
  int c=memcmp(&t1,&t2,SB);*t=(c<=0)?t1:t2;}
St*states;size_t nst=0,cap;uint32_t*htab;size_t hcap;
static size_t find_or_add(const St*s){uint64_t h=hsh(s);size_t i=h&(hcap-1);
  while(htab[i]!=0xffffffffu){if(eq(&states[htab[i]],s))return htab[i];i=(i+1)&(hcap-1);}
  if(nst>=cap){cap=cap*5/4+16;states=realloc(states,cap*sizeof(St));}
  states[nst]=*s;htab[i]=nst;return nst++;}
uint64_t*eoff;uint32_t*edst;uint8_t*etag;size_t nedges=0,ecap;
int main(int argc,char**argv){
  K=atoi(argv[1]);int iters=atoi(argv[2]);TOT=atoi(argv[3]);int lh=atoi(argv[4]);size_t scap=(size_t)atoll(argv[5]);
  cap=scap;states=malloc(cap*sizeof(St));hcap=1ULL<<lh;htab=malloc(hcap*4);memset(htab,0xff,hcap*4);
  St init;memset(&init,0,sizeof init);find_or_add(&init);
  for(int pass=0;pass<2;pass++){
  if(pass==1){eoff=malloc((nst+1)*8);ecap=nst*4;edst=malloc(ecap*4);etag=malloc(ecap);nedges=0;fprintf(stderr,"states %zu; allocating %zu edge slots\n",nst,ecap);}
  for(size_t i=0;i<nst;i++){
    if(pass==1)eoff[i]=nedges;St s=states[i];int r=s.r,l=s.l;
    for(int er=0;er<2;er++)for(int el=0;el<2;el++){
      uint8_t R[MAXW+2],L[MAXW+2];int nr=r,nl=l;for(int q=0;q<r;q++)R[q]=getl(&s,q);for(int q=0;q<l;q++)L[q]=getl(&s,r+q);
      int closed[2],nc=0; // label or 0 for forgotten
      if(!er){closed[nc++]=nr?R[--nr]:0;}
      if(!el){closed[nc++]=nl?L[--nl]:0;}
      if(nc==2&&closed[0]&&closed[0]==closed[1])continue; // cycle among known arcs
      int cur=0;
      if(nc==2){int a=closed[0],b=closed[1];if(a&&b){for(int q=0;q<nr;q++)if(R[q]==b)R[q]=a;for(int q=0;q<nl;q++)if(L[q]==b)L[q]=a;cur=a;}}
      else if(nc==1)cur=closed[0];
      if(er||el){if(!cur){int mx=0;for(int q=0;q<nr;q++)if(R[q]>mx)mx=R[q];for(int q=0;q<nl;q++)if(L[q]>mx)mx=L[q];cur=mx+1;}
        if(er)R[nr++]=cur;if(el)L[nl++]=cur;}
      // singleton labels carry no information: replace by '?' (0); strip '?' from the bottom of each stack
      {int cnt[64];memset(cnt,0,sizeof cnt);for(int q=0;q<nr;q++)cnt[R[q]]++;for(int q=0;q<nl;q++)cnt[L[q]]++;
       for(int q=0;q<nr;q++)if(R[q]&&cnt[R[q]]==1)R[q]=0;for(int q=0;q<nl;q++)if(L[q]&&cnt[L[q]]==1)L[q]=0;
       int a=0;while(a<nr&&R[a]==0)a++;memmove(R,R+a,nr-a);nr-=a;a=0;while(a<nl&&L[a]==0)a++;memmove(L,L+a,nl-a);nl-=a;}
      // forgetting
      if(TOT>0){ // total budget: drop deepest arcs, alternating sides, from the taller side first
        while(nr+nl>TOT){if(nr>=nl){memmove(R,R+1,--nr);}else{memmove(L,L+1,--nl);}}
      }else{while(nr>K){memmove(R,R+1,--nr);}while(nl>K){memmove(L,L+1,--nl);}}
      for(int rep=0;rep<2;rep++){int cnt[64];memset(cnt,0,sizeof cnt);for(int q=0;q<nr;q++)cnt[R[q]]++;for(int q=0;q<nl;q++)cnt[L[q]]++;
       for(int q=0;q<nr;q++)if(R[q]&&cnt[R[q]]==1)R[q]=0;for(int q=0;q<nl;q++)if(L[q]&&cnt[L[q]]==1)L[q]=0;
       int a=0;while(a<nr&&R[a]==0)a++;memmove(R,R+a,nr-a);nr-=a;a=0;while(a<nl&&L[a]==0)a++;memmove(L,L+a,nl-a);nl-=a;}
      St t;make_state(&t,R,nr,L,nl);
      if(nst*3>hcap*2){fprintf(stderr,"hash too small\n");return 1;}
      if(nst>=cap){fprintf(stderr,"STATECAP too small\n");return 1;}
      size_t j=find_or_add(&t);
      if(pass==1){etag[nedges]=((nc==2&&(closed[0]==0||closed[1]==0))||(!er&&r==0)||(!el&&l==0))?1:0;edst[nedges++]=(uint32_t)j;}}}
  }
  eoff[nst]=nedges;free(htab);free(states);
  fprintf(stderr,"k=%d tot=%d states=%zu edges=%zu\n",K,TOT,nst,nedges);
  size_t N=nst;double*x=malloc(N*8),*y=malloc(N*8);for(size_t i=0;i<N;i++)x[i]=1;double l2=0;
  for(int it=0;it<iters;it++){for(int h=0;h<2;h++){
    #pragma omp parallel for schedule(static)
    for(size_t i=0;i<N;i++){double s=0;for(uint64_t k=eoff[i];k<eoff[i+1];k++)if(!etag[k])s+=x[edst[k]];y[i]=s;}
    double mx=0;for(size_t i=0;i<N;i++)if(y[i]>mx)mx=y[i];
    #pragma omp parallel for schedule(static)
    for(size_t i=0;i<N;i++)x[i]=y[i]/mx;if(h==0)l2=mx;else l2*=mx;}
    if(it%50==49)fprintf(stderr,"it %d lam_restricted %.9f\n",it+1,l2);}
  int64_t*w=(int64_t*)x,*t1=malloc(N*8),*t2=(int64_t*)y;
  for(size_t i=0;i<N;i++){int64_t v=(int64_t)floor(x[i]*(double)(1LL<<40));w[i]=v<0?0:v;}
  #pragma omp parallel for schedule(static)
  for(size_t i=0;i<N;i++){int64_t s=0;for(uint64_t k=eoff[i];k<eoff[i+1];k++)if(!etag[k])s+=w[edst[k]];t1[i]=s;}
  #pragma omp parallel for schedule(static)
  for(size_t i=0;i<N;i++){int64_t s=0;for(uint64_t k=eoff[i];k<eoff[i+1];k++)if(!etag[k])s+=t1[edst[k]];t2[i]=s;}
  uint8_t*seen=calloc(N,1);uint32_t*q=malloc(N*4);size_t qh=0,qt=0;seen[0]=1;q[qt++]=0;
  while(qh<qt){uint32_t i=q[qh++];for(uint64_t k=eoff[i];k<eoff[i+1];k++)if(!etag[k]){uint32_t j=edst[k];if(!seen[j]){seen[j]=1;q[qt++]=j;}}}
  int64_t bn=0,bd=1;int first=1;size_t npos=0;
  for(size_t i=0;i<N;i++){if(w[i]<=0||!seen[i])continue;npos++;if(first){bn=t2[i];bd=w[i];first=0;}else if((__int128)t2[i]*bd<(__int128)bn*w[i]){bn=t2[i];bd=w[i];}}
  printf("RESTRICTED tot=%d states=%zu reachable=%zu rho_est(T''^2)=%.6f CERTIFIED gamma^2 >= %lld/%lld = %.6f (positive+reachable %zu)\n",TOT,nst,qt,l2,(long long)bn,(long long)bd,(double)bn/(double)bd,npos);
  return 0;}
