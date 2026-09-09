// Upper bound on the meander constant R from the label-forgetting relaxation (per-side memory k).
// State: upper/lower stacks of labelled arcs (at most k per side; deeper arcs forgotten). Every closed meander
// is a walk of this automaton, so rho(T)^2 >= R.  Certificate: integer u>0 with T^2 u <= c u entrywise => rho(T^2) <= c.
// Usage: relaxc k iters [total]   (total>0: keep at most `total` labelled arcs overall instead of k per side)
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#define MAXW 28
typedef struct { uint8_t r,l; uint8_t p[14]; } St; // 16 bytes, nibble-packed labels 0..15 (0 = forgotten)
int K,TOT;
static uint64_t hsh(const St*a){uint64_t h=1469598103934665603ULL;const uint8_t*b=(const uint8_t*)a;for(int i=0;i<16;i++){h^=b[i];h*=1099511628211ULL;}h^=h>>29;h*=0x9E3779B97F4A7C15ULL;h^=h>>32;return h;}
static int eq(const St*a,const St*b){return memcmp(a,b,16)==0;}
static inline int getl(const St*s,int i){return (s->p[i>>1]>>((i&1)*4))&15;}
static void pack(St*s,const uint8_t*lab,int n){memset(s->p,0,14);for(int i=0;i<n;i++)s->p[i>>1]|=lab[i]<<((i&1)*4);}
static void canon_lab(uint8_t*lab,int n){uint8_t map[64];memset(map,0,64);int nxt=1;for(int i=0;i<n;i++){uint8_t c=lab[i];if(!c)continue;if(!map[c])map[c]=nxt++;lab[i]=map[c];}}
static void make_state(St*t,const uint8_t*R,int nr,const uint8_t*L,int nl){
  uint8_t a[MAXW+2],b[MAXW+2];memcpy(a,R,nr);memcpy(a+nr,L,nl);canon_lab(a,nr+nl);memcpy(b,L,nl);memcpy(b+nl,R,nr);canon_lab(b,nr+nl);
  St t1,t2;memset(&t1,0,16);memset(&t2,0,16);t1.r=nr;t1.l=nl;pack(&t1,a,nr+nl);t2.r=nl;t2.l=nr;pack(&t2,b,nr+nl);
  int c=memcmp(&t1,&t2,16);*t=(c<=0)?t1:t2;}
St*states;size_t nst=0,cap;uint32_t*htab;size_t hcap;
static size_t find_or_add(const St*s){uint64_t h=hsh(s);size_t i=h&(hcap-1);
  while(htab[i]!=0xffffffffu){if(eq(&states[htab[i]],s))return htab[i];i=(i+1)&(hcap-1);}
  if(nst>=cap){cap=cap*5/4+16;states=realloc(states,cap*sizeof(St));}
  states[nst]=*s;htab[i]=nst;return nst++;}
uint32_t*eoff,*edst;size_t nedges=0,ecap;
int main(int argc,char**argv){
  K=atoi(argv[1]);int iters=atoi(argv[2]);TOT=argc>3?atoi(argv[3]):0;
  cap=1<<20;states=malloc(cap*sizeof(St));hcap=1ULL<<28;htab=malloc(hcap*4);memset(htab,0xff,hcap*4);
  ecap=1<<22;edst=malloc(ecap*4);
  St init;memset(&init,0,sizeof init);find_or_add(&init);
  size_t eoffcap=cap;eoff=malloc((eoffcap+1)*4);
  for(size_t i=0;i<nst;i++){
    if(i+1>=eoffcap){eoffcap=eoffcap*5/4;eoff=realloc(eoff,(eoffcap+1)*4);}
    eoff[i]=nedges;St s=states[i];int r=s.r,l=s.l;
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
      if(nst*2>hcap){fprintf(stderr,"hash too small\n");return 1;}
      size_t j=find_or_add(&t);
      if(nedges>=ecap){ecap=ecap*5/4;edst=realloc(edst,ecap*4);}
      edst[nedges++]=j;}}
  eoff[nst]=nedges;free(htab);free(states);
  fprintf(stderr,"k=%d tot=%d states=%zu edges=%zu\n",K,TOT,nst,nedges);
  double*v=malloc(nst*8),*u=malloc(nst*8);for(size_t i=0;i<nst;i++)v[i]=1;
  double lam=0;
  for(int it=0;it<iters;it++){ // right vector: v <- T^2 v  (v_i = sum_j T_ij v_j)
    for(int h=0;h<2;h++){for(size_t i=0;i<nst;i++){double s=0;for(uint32_t k=eoff[i];k<eoff[i+1];k++)s+=v[edst[k]];u[i]=s;}
      double mx=0;for(size_t i=0;i<nst;i++)if(u[i]>mx)mx=u[i];for(size_t i=0;i<nst;i++)v[i]=u[i]/mx;if(h==0)lam=mx;else lam*=mx;}
    if(it%50==49)fprintf(stderr,"it %d lam %.9f\n",it+1,lam);}
  // certificate: integer w = ceil(v*2^40) with floor 1; c = max_i (T^2 w)_i / w_i
  int64_t*w=malloc(nst*8),*t1=(int64_t*)u,*t2;
  for(size_t i=0;i<nst;i++){int64_t x=(int64_t)ceil(v[i]*(double)(1LL<<40));if(x<1)x=1;w[i]=x;}
  t2=(int64_t*)v;
  for(size_t i=0;i<nst;i++){int64_t s=0;for(uint32_t k=eoff[i];k<eoff[i+1];k++)s+=w[edst[k]];t1[i]=s;}
  for(size_t i=0;i<nst;i++){int64_t s=0;for(uint32_t k=eoff[i];k<eoff[i+1];k++)s+=t1[edst[k]];t2[i]=s;}
  int64_t bn=0,bd=1;for(size_t i=0;i<nst;i++){if((__int128)t2[i]*bd>(__int128)bn*w[i]){bn=t2[i];bd=w[i];}}
  double c=(double)bn/(double)bd;
  printf("k=%d tot=%d states=%zu rho_est(T^2)=%.6f CERTIFIED rho(T^2) <= %lld/%lld = %.6f => R <= %.5f\n",K,TOT,nst,lam,(long long)bn,(long long)bd,c,c*(1+1e-15));
  return 0;}
