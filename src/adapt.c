// Adaptive (Perron-mass-selected) truncation of the closed-meander transfer matrix, mirror quotient.
// Usage: adapt W iters K1 K2 ...   (K = number of orbit-states to keep; K<0 means "width <= -K" baseline)
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
int W;
typedef struct { uint8_t r,l; uint8_t p[14]; } St;
static inline int getl(const St*s,int i){return (s->p[i>>1]>>((i&1)*4))&15;}
static void pack(St*s,const uint8_t*lab,int n){memset(s->p,0,14);for(int i=0;i<n;i++)s->p[i>>1]|=lab[i]<<((i&1)*4);}
static uint64_t hsh(const St*a){uint64_t h=1469598103934665603ULL;const uint8_t*b=(const uint8_t*)a;for(int i=0;i<16;i++){h^=b[i];h*=1099511628211ULL;}return h;}
St*states;size_t nst=0,cap;uint32_t*htab;size_t hcap;
size_t find_or_add(const St*s){uint64_t h=hsh(s);size_t i=h&(hcap-1);
  while(htab[i]!=0xffffffffu){if(memcmp(&states[htab[i]],s,16)==0)return htab[i];i=(i+1)&(hcap-1);}
  if(nst>=cap){cap=cap*5/4;states=realloc(states,cap*sizeof(St));}
  states[nst]=*s;htab[i]=nst;return nst++;}
static void canon_lab(uint8_t*lab,int n){uint8_t map[32];memset(map,0,32);int nxt=1;for(int i=0;i<n;i++){uint8_t c=lab[i];if(!map[c])map[c]=nxt++;lab[i]=map[c];}}
static void make_state(St*t,uint8_t*R,int nr,uint8_t*L,int nl){
  uint8_t a[32],b[32];memcpy(a,R,nr);memcpy(a+nr,L,nl);canon_lab(a,nr+nl);
  memcpy(b,L,nl);memcpy(b+nl,R,nr);canon_lab(b,nr+nl);
  St t1,t2;t1.r=nr;t1.l=nl;pack(&t1,a,nr+nl);t2.r=nl;t2.l=nr;pack(&t2,b,nr+nl);
  *t=(memcmp(&t1,&t2,16)<=0)?t1:t2;}
uint32_t*eoff,*edst;size_t nedges=0,ecap;
uint8_t*width;

static void build(void){
  St init;memset(&init,0,16);find_or_add(&init);
  size_t eoffcap=cap;eoff=malloc(eoffcap*4);
  for(size_t i=0;i<nst;i++){
    if(i+1>=eoffcap){eoffcap=eoffcap*5/4;eoff=realloc(eoff,eoffcap*4);}
    eoff[i]=nedges;St s=states[i];int r=s.r,l=s.l;
    uint8_t R0[32],L0[32];for(int k=0;k<r;k++)R0[k]=getl(&s,k);for(int k=0;k<l;k++)L0[k]=getl(&s,r+k);
    for(int er=0;er<2;er++)for(int el=0;el<2;el++){
      if(!er&&r==0)continue;if(!el&&l==0)continue;
      if(r+l+(er?1:-1)+(el?1:-1)>W)continue;
      uint8_t R[32],L[32];int nr=r,nl=l;memcpy(R,R0,r);memcpy(L,L0,l);
      int touched[2],nt=0;if(!er)touched[nt++]=R[--nr];if(!el)touched[nt++]=L[--nl];
      if(nt==2&&touched[0]==touched[1])continue;
      int cid;
      if(nt==0){int mx=0;for(int k=0;k<nr;k++)if(R[k]>mx)mx=R[k];for(int k=0;k<nl;k++)if(L[k]>mx)mx=L[k];cid=mx+1;}
      else{cid=touched[0];if(nt==2&&touched[1]<cid)cid=touched[1];
        for(int k=0;k<nr;k++)for(int q=0;q<nt;q++)if(R[k]==touched[q])R[k]=cid;
        for(int k=0;k<nl;k++)for(int q=0;q<nt;q++)if(L[k]==touched[q])L[k]=cid;}
      if(er)R[nr++]=cid;if(el)L[nl++]=cid;
      St t;make_state(&t,R,nr,L,nl);
      if(nst*2>hcap){fprintf(stderr,"hash too small\n");exit(1);}
      size_t j=find_or_add(&t);
      if(nedges>=ecap){ecap=ecap*5/4;edst=realloc(edst,ecap*4);}
      edst[nedges++]=j;}}
  eoff[nst]=nedges;
  width=malloc(nst);for(size_t i=0;i<nst;i++)width[i]=states[i].r+states[i].l;
}

// Perron vectors of T restricted to the subset flagged by keep[] (keep==NULL: all). Returns rho estimate of T^2.
// a: left (a^T T = rho a^T, i.e. push forward), b: right (T b = rho b, pull back).
static double perron(const uint8_t*keep,int iters,double*a,double*b,double*u){
  size_t n=nst;
  for(size_t i=0;i<n;i++){a[i]=(keep&&!keep[i])?0:1;b[i]=a[i];}
  double lam=1;
  for(int it=0;it<iters;it++){
    for(int half=0;half<2;half++){
      memset(u,0,n*8);
      for(size_t i=0;i<n;i++){double vi=a[i];if(vi==0)continue;for(uint32_t k=eoff[i];k<eoff[i+1];k++){uint32_t j=edst[k];if(!keep||keep[j])u[j]+=vi;}}
      double mx=0;for(size_t i=0;i<n;i++)if(u[i]>mx)mx=u[i];
      for(size_t i=0;i<n;i++)a[i]=u[i]/mx;
      if(half==0)lam=mx;
      else lam*=mx;
    }
    for(int half=0;half<2;half++){
      memset(u,0,n*8);
      for(size_t i=0;i<n;i++){if(keep&&!keep[i])continue;double s=0;for(uint32_t k=eoff[i];k<eoff[i+1];k++){uint32_t j=edst[k];s+=b[j];}u[i]=s;}
      double mx=0;for(size_t i=0;i<n;i++)if(u[i]>mx)mx=u[i];
      for(size_t i=0;i<n;i++)b[i]=u[i]/mx;
    }
  }
  return lam;
}

// exact integer certificate on subset: min over positive entries of ((T^T)^2 v)_i / v_i, using left vector a
static double certify(const uint8_t*keep,const double*a,int64_t*vi,int64_t*t1,int64_t*t2,long long*pn,long long*pd){
  size_t n=nst;
  for(size_t i=0;i<n;i++){double q=(keep&&!keep[i])?0:a[i]*(double)(1LL<<52);vi[i]=(int64_t)floor(q);if(vi[i]<0)vi[i]=0;}
  memset(t1,0,n*8);for(size_t i=0;i<n;i++){int64_t x=vi[i];if(x)for(uint32_t k=eoff[i];k<eoff[i+1];k++){uint32_t j=edst[k];if(!keep||keep[j])t1[j]+=x;}}
  memset(t2,0,n*8);for(size_t i=0;i<n;i++){int64_t x=t1[i];if(x)for(uint32_t k=eoff[i];k<eoff[i+1];k++){uint32_t j=edst[k];if(!keep||keep[j])t2[j]+=x;}}
  int64_t bn=0,bd=1;int first=1;
  for(size_t i=0;i<n;i++){if(vi[i]<=0)continue;
    if(first){bn=t2[i];bd=vi[i];first=0;}else if((__int128)t2[i]*bd<(__int128)bn*vi[i]){bn=t2[i];bd=vi[i];}}
  *pn=bn;*pd=bd;return (double)bn/(double)bd;
}

static double*gmass;
static int cmpm(const void*x,const void*y){double a=gmass[*(const uint32_t*)x],b=gmass[*(const uint32_t*)y];return a<b?1:(a>b?-1:0);}

int main(int argc,char**argv){
  W=atoi(argv[1]);int iters=atoi(argv[2]);
  cap=1<<20;states=malloc(cap*sizeof(St));hcap=1ULL<<26;htab=malloc(hcap*4);memset(htab,0xff,hcap*4);
  ecap=1<<22;edst=malloc(ecap*4);
  build();free(htab);
  fprintf(stderr,"W=%d orbit-states=%zu edges=%zu\n",W,nst,nedges);
  double*a=malloc(nst*8),*b=malloc(nst*8),*u=malloc(nst*8);
  double rho=perron(NULL,iters,a,b,u);
  fprintf(stderr,"full: rho(T^2)~%.6f\n",rho);
  // Perron vectors of T (not T^2): a + T^T a / sqrt(rho) etc.
  double sr=sqrt(rho);
  {memset(u,0,nst*8);for(size_t i=0;i<nst;i++){double vi=a[i];for(uint32_t k=eoff[i];k<eoff[i+1];k++)u[edst[k]]+=vi;}for(size_t i=0;i<nst;i++)a[i]+=u[i]/sr;}
  {memset(u,0,nst*8);for(size_t i=0;i<nst;i++){double s=0;for(uint32_t k=eoff[i];k<eoff[i+1];k++)s+=b[edst[k]];u[i]=s;}for(size_t i=0;i<nst;i++)b[i]+=u[i]/sr;}
  double*mass=malloc(nst*8);double tot=0;for(size_t i=0;i<nst;i++){mass[i]=a[i]*b[i];tot+=mass[i];}
  double byw[64]={0};size_t cw[64]={0};for(size_t i=0;i<nst;i++){byw[width[i]]+=mass[i];cw[width[i]]++;}
  fprintf(stderr,"mass by width:");for(int w=0;w<=W;w++)if(cw[w])fprintf(stderr," %d:%.2f%%(%zu)",w,100*byw[w]/tot,cw[w]);fprintf(stderr,"\n");
  uint32_t*ord=malloc(nst*4);for(size_t i=0;i<nst;i++)ord[i]=i;gmass=mass;qsort(ord,nst,4,cmpm);
  uint8_t*keep=malloc(nst);
  double*a2=malloc(nst*8),*b2=malloc(nst*8);
  int64_t*vi=malloc(nst*8),*t1=malloc(nst*8),*t2=malloc(nst*8);
  for(int ai=3;ai<argc;ai++){
    long K=atol(argv[ai]);size_t cnt=0;
    if(K<0){int Wp=-K;for(size_t i=0;i<nst;i++){keep[i]=width[i]<=Wp;cnt+=keep[i];}}
    else{memset(keep,0,nst);for(long i=0;i<K&&i<(long)nst;i++){keep[ord[i]]=1;cnt++;}}
    // ensure connectivity from init(0)? init has width 0 and is in every width set; in mass sets it may be absent -- fine (dropped in muq too).
    double r2=perron(keep,iters,a2,b2,u);
    long long pn,pd;double c=certify(keep,a2,vi,t1,t2,&pn,&pd);
    double maxw=0;size_t nw[64]={0};for(size_t i=0;i<nst;i++)if(keep[i])nw[width[i]]++;
    printf("K=%ld states=%zu rho_est=%.6f certified=%lld/%lld=%.6f  widths:",K,cnt,r2,pn,pd,c);
    for(int w=0;w<=W;w++)if(nw[w])printf(" %d:%zu",w,nw[w]);printf("\n");fflush(stdout);
  }
  return 0;}
