// Adaptive ("selected-configuration") truncation of the closed-meander transfer matrix, mirror quotient.
// Start: all orbit-states of width <= W0.  Then repeatedly: Perron vectors on S -> score frontier states
// (t not in S) by  a_est(t)*b_est(t), a_est = sum_{i in S} a_i T_it / rho, b_est = sum_{j in S} T_tj b_j / rho,
// keep the top-K of S u F, rebuild, iterate.  At each stage an exact integer certificate
// (T_S^T)^2 v >= c v is computed; every such c is a rigorous lower bound on R.
// Usage: grow W0 WMAX ITERS OUTPREFIX K1 K2 K3 ...
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <time.h>
#define SB 20
typedef struct { uint8_t r,l; uint8_t p[SB-2]; } St; // labels 1..16 stored as nibbles 0..15; width <= 36 in bytes, <=32 paths-wise
static inline int getl(const St*s,int i){return ((s->p[i>>1]>>((i&1)*4))&15)+1;}
static void pack(St*s,const uint8_t*lab,int n){memset(s->p,0,SB-2);for(int i=0;i<n;i++)s->p[i>>1]|=(lab[i]-1)<<((i&1)*4);}
static inline uint64_t hsh(const St*a){uint64_t h=1469598103934665603ULL;const uint8_t*b=(const uint8_t*)a;for(int i=0;i<SB;i++){h^=b[i];h*=1099511628211ULL;}h^=h>>29;h*=0x9E3779B97F4A7C15ULL;h^=h>>32;return h;}
static void canon_lab(uint8_t*lab,int n){uint8_t map[40];memset(map,0,40);int nxt=1;for(int i=0;i<n;i++){uint8_t c=lab[i];if(!map[c])map[c]=nxt++;lab[i]=map[c];}}
static void make_state(St*t,uint8_t*R,int nr,uint8_t*L,int nl){
  uint8_t a[40],b[40];memcpy(a,R,nr);memcpy(a+nr,L,nl);canon_lab(a,nr+nl);
  memcpy(b,L,nl);memcpy(b+nl,R,nr);canon_lab(b,nr+nl);
  St t1,t2;memset(&t1,0,SB);memset(&t2,0,SB);t1.r=nr;t1.l=nl;pack(&t1,a,nr+nl);t2.r=nl;t2.l=nr;pack(&t2,b,nr+nl);
  *t=(memcmp(&t1,&t2,SB)<=0)?t1:t2;}
static double now(void){struct timespec t;clock_gettime(CLOCK_MONOTONIC,&t);return t.tv_sec+1e-9*t.tv_nsec;}

int WMAX;
// successors of s with width <= wcap; returns count, fills out[]
static int succ(const St*s,int wcap,St*out){
  int r=s->r,l=s->l,n=0;uint8_t R0[40],L0[40];for(int k=0;k<r;k++)R0[k]=getl(s,k);for(int k=0;k<l;k++)L0[k]=getl(s,r+k);
  for(int er=0;er<2;er++)for(int el=0;el<2;el++){
    if(!er&&r==0)continue;if(!el&&l==0)continue;
    if(r+l+(er?1:-1)+(el?1:-1)>wcap)continue;
    uint8_t R[40],L[40];int nr=r,nl=l;memcpy(R,R0,r);memcpy(L,L0,l);
    int touched[2],nt=0;if(!er)touched[nt++]=R[--nr];if(!el)touched[nt++]=L[--nl];
    if(nt==2&&touched[0]==touched[1])continue;
    int cid;
    if(nt==0){int mx=0;for(int k=0;k<nr;k++)if(R[k]>mx)mx=R[k];for(int k=0;k<nl;k++)if(L[k]>mx)mx=L[k];cid=mx+1;}
    else{cid=touched[0];if(nt==2&&touched[1]<cid)cid=touched[1];
      for(int k=0;k<nr;k++)for(int q=0;q<nt;q++)if(R[k]==touched[q])R[k]=cid;
      for(int k=0;k<nl;k++)for(int q=0;q<nt;q++)if(L[k]==touched[q])L[k]=cid;}
    if(er)R[nr++]=cid;if(el)L[nl++]=cid;
    make_state(&out[n++],R,nr,L,nl);}
  return n;}

// ---- the current set S ----
St*states;size_t nst=0,cap;uint32_t*htab;size_t hcap;
static void hash_clear(void){memset(htab,0xff,hcap*4);}
static void hash_insert_all(void){hash_clear();for(size_t i=0;i<nst;i++){uint64_t h=hsh(&states[i]);size_t k=h&(hcap-1);while(htab[k]!=0xffffffffu)k=(k+1)&(hcap-1);htab[k]=i;}}
static inline int64_t lookup(const St*s){uint64_t h=hsh(s);size_t k=h&(hcap-1);while(htab[k]!=0xffffffffu){if(memcmp(&states[htab[k]],s,SB)==0)return htab[k];k=(k+1)&(hcap-1);}return -1;}
static size_t find_or_add(const St*s){uint64_t h=hsh(s);size_t k=h&(hcap-1);
  while(htab[k]!=0xffffffffu){if(memcmp(&states[htab[k]],s,SB)==0)return htab[k];k=(k+1)&(hcap-1);}
  if(nst>=cap){cap=cap*5/4+16;states=realloc(states,cap*sizeof(St));}
  states[nst]=*s;htab[k]=nst;return nst++;}
uint32_t*eoff=NULL,*edst=NULL;size_t nedges=0,ecap=0;

static void build_csr(void){ // transitions inside S (width cap WMAX)
  if(!eoff)eoff=malloc((cap+1)*4);else eoff=realloc(eoff,(cap+1)*4);
  nedges=0;St out[4];
  for(size_t i=0;i<nst;i++){eoff[i]=nedges;int n=succ(&states[i],WMAX,out);
    for(int q=0;q<n;q++){int64_t j=lookup(&out[q]);if(j<0)continue;
      if(nedges>=ecap){ecap=ecap*5/4+1024;edst=realloc(edst,ecap*4);}edst[nedges++]=(uint32_t)j;}}
  eoff[nst]=nedges;}

// power iteration on T_S^2 for left vector a (warm start), plus right vector b. u scratch. returns rho(T^2) est.
static double perron(int iters,float*a,float*b,float*u,double*lastchange){
  size_t n=nst;double lam=0,l1=0,l2=0;
  for(int it=0;it<iters;it++){
    for(int half=0;half<2;half++){
      memset(u,0,n*4);
      for(size_t i=0;i<n;i++){float vi=a[i];if(vi==0)continue;for(uint32_t k=eoff[i];k<eoff[i+1];k++)u[edst[k]]+=vi;}
      double mx=0;for(size_t i=0;i<n;i++)if(u[i]>mx)mx=u[i];
      for(size_t i=0;i<n;i++)a[i]=u[i]/mx;
      if(half==0)l1=mx;else l2=mx;}
    for(int half=0;half<2;half++){
      memset(u,0,n*4);
      for(size_t i=0;i<n;i++){double s=0;for(uint32_t k=eoff[i];k<eoff[i+1];k++)s+=b[edst[k]];u[i]=s;}
      double mx=0;for(size_t i=0;i<n;i++)if(u[i]>mx)mx=u[i];
      for(size_t i=0;i<n;i++)b[i]=u[i]/mx;}
    double nl=l1*l2;if(lastchange)*lastchange=fabs(nl-lam);lam=nl;
    if(it%25==24)fprintf(stderr,"   it %d lam %.8f\n",it+1,lam);}
  return lam;}

// exact certificate: v = floor(a*2^52) (double from float), t2=(T^T)^2 v, c = min t2_i/v_i over v_i>0
static double certify(const float*a,long long*pn,long long*pd,int64_t**pv){
  size_t n=nst;int64_t*vi=malloc(n*8),*t1=calloc(n,8),*t2=calloc(n,8);
  for(size_t i=0;i<n;i++){double q=(double)a[i]*(double)(1LL<<52);vi[i]=(int64_t)floor(q);if(vi[i]<0)vi[i]=0;}
  for(size_t i=0;i<n;i++){int64_t x=vi[i];if(x)for(uint32_t k=eoff[i];k<eoff[i+1];k++)t1[edst[k]]+=x;}
  for(size_t i=0;i<n;i++){int64_t x=t1[i];if(x)for(uint32_t k=eoff[i];k<eoff[i+1];k++)t2[edst[k]]+=x;}
  int64_t bn=0,bd=1;int first=1;size_t npos=0;
  for(size_t i=0;i<n;i++){if(vi[i]<=0)continue;npos++;
    if(first){bn=t2[i];bd=vi[i];first=0;}else if((__int128)t2[i]*bd<(__int128)bn*vi[i]){bn=t2[i];bd=vi[i];}}
  free(t1);free(t2);*pn=bn;*pd=bd;if(pv)*pv=vi;else free(vi);
  fprintf(stderr,"   certificate positive entries %zu of %zu\n",npos,n);
  return (double)bn/(double)bd;}

// ---- frontier hash (chunked) ----
typedef struct { St s; float aest; } FE; // 24 bytes
FE*ftab;size_t fcap;
static inline void fclear(void){for(size_t i=0;i<fcap;i++)ftab[i].s.r=255;}
static size_t fcount;
static inline FE* fget(const St*s){uint64_t h=hsh(s);size_t k=(h>>20)&(fcap-1);while(ftab[k].s.r!=255){if(memcmp(&ftab[k].s,s,SB)==0)return &ftab[k];k=(k+1)&(fcap-1);}
  if(fcount*10>fcap*8){fprintf(stderr,"frontier hash full\n");exit(1);}ftab[k].s=*s;ftab[k].aest=0;fcount++;return &ftab[k];}

static int cmpf(const void*x,const void*y){float a=*(const float*)x,b=*(const float*)y;return a<b?1:(a>b?-1:0);}

int main(int argc,char**argv){
  int W0=atoi(argv[1]);WMAX=atoi(argv[2]);int iters=atoi(argv[3]);const char*outp=argv[4];
  cap=1<<22;states=malloc(cap*sizeof(St));hcap=1ULL<<27;htab=malloc(hcap*4);hash_clear();
  fcap=1ULL<<25;ftab=malloc(fcap*sizeof(FE));
  // initial S: BFS width<=W0
  float*a,*b,*u;
  if(W0>=0){St init;memset(&init,0,SB);find_or_add(&init);St out[4];
   for(size_t i=0;i<nst;i++){int n=succ(&states[i],W0,out);for(int q=0;q<n;q++)find_or_add(&out[q]);}
   fprintf(stderr,"initial S: width<=%d, %zu orbit-states\n",W0,nst);
   build_csr();a=malloc(nst*4);b=malloc(nst*4);for(size_t i=0;i<nst;i++){a[i]=1;b[i]=1;}}
  else{char fn[512];snprintf(fn,512,"%s_ck_states.bin",outp);FILE*f=fopen(fn,"rb");fseek(f,0,SEEK_END);nst=ftell(f)/SB;fseek(f,0,SEEK_SET);
   cap=nst+1024;states=realloc(states,cap*SB);if(fread(states,SB,nst,f)!=nst)return 1;fclose(f);
   a=malloc(nst*4);b=malloc(nst*4);snprintf(fn,512,"%s_ck_ab.bin",outp);f=fopen(fn,"rb");if(fread(a,4,nst,f)!=nst||fread(b,4,nst,f)!=nst)return 1;fclose(f);
   hash_insert_all();build_csr();fprintf(stderr,"resumed %zu states\n",nst);}
  u=malloc(nst*4);
  St j0;{uint8_t R[1]={1},L[1]={1};make_state(&j0,R,1,L,1);}
  for(int stage=0;stage+5<=argc;stage++){
    double t0=now();
    int itn=(stage+5>=argc||atol(argv[stage+5])>0)?iters:iters/3+1;
    double rho=perron(itn,a,b,u,NULL);double sr=sqrt(rho);
    long long pn,pd;double c=certify(a,&pn,&pd,NULL);
    size_t nw[64]={0};for(size_t i=0;i<nst;i++)nw[states[i].r+states[i].l]++;
    printf("STAGE %d: states=%zu edges=%zu rho_est=%.8f CERTIFIED rho(T^2) >= %lld/%lld = %.8f  [%.0fs]\n  widths:",stage,nst,nedges,rho,pn,pd,c,now()-t0);
    for(int w=0;w<64;w++)if(nw[w])printf(" %d:%zu",w,nw[w]);printf("\n");fflush(stdout);
    {int64_t jj=lookup(&j0);if(jj<0||a[jj]<=0)printf("  WARNING: j0 missing or zero\n");}
    {char fn[512];snprintf(fn,512,"%s_ck_states.bin",outp);FILE*f=fopen(fn,"wb");fwrite(states,SB,nst,f);fclose(f);
     snprintf(fn,512,"%s_ck_ab.bin",outp);f=fopen(fn,"wb");fwrite(a,4,nst,f);fwrite(b,4,nst,f);fclose(f);}
    if(stage+5>=argc)break;
    long K=atol(argv[stage+5]);
    // T-Perron vectors: aT = a + T^T a / sr ; bT = b + T b / sr (store in a,b temporarily via u)
    memset(u,0,nst*4);for(size_t i=0;i<nst;i++){float vi=a[i];for(uint32_t k=eoff[i];k<eoff[i+1];k++)u[edst[k]]+=vi;}
    float*aT=malloc(nst*4),*bT=malloc(nst*4);for(size_t i=0;i<nst;i++)aT[i]=a[i]+u[i]/sr;
    for(size_t i=0;i<nst;i++){double s=0;for(uint32_t k=eoff[i];k<eoff[i+1];k++)s+=b[edst[k]];bT[i]=b[i]+s/sr;}
    // scores of S members
    float*sc=malloc(nst*4);for(size_t i=0;i<nst;i++)sc[i]=aT[i]*bT[i];
    // pass 1: frontier scores (chunked)
    int P=(int)(3.0*nst/(0.75*fcap))+1;
    #define NB 16384
    static size_t hist[NB];memset(hist,0,sizeof hist);size_t nfs=0;
    #define BIN(x) ({int _b=(int)floor(64.0*log2((double)(x)))+NB-64; _b<0?0:(_b>=NB?NB-1:_b);})
    for(size_t i=0;i<nst;i++)if(sc[i]>0)hist[BIN(sc[i])]++;
    St out[4],out2[4];
    for(int pass=0;pass<2;pass++){
      float tau=0;
      if(pass==1){ // threshold: K-th largest of all scores
        if(K<0)tau=-1;else{size_t cum=0;int bb=NB-1;for(;bb>0;bb--){cum+=hist[bb];if(cum>=(size_t)K)break;}
          tau=(float)exp2((bb-NB+64)/64.0);if(cum<(size_t)K)tau=0;}
        fprintf(stderr,"  frontier %zu, tau=%g\n",nfs,tau);}
      // pass 1 collects fs; pass 2 collects new states into newS
      St*newS=NULL;size_t nnew=0,ncap=0;
      if(pass==1){ncap=(K<0?nst+nfs:K)+1024;newS=malloc(ncap*sizeof(St));for(size_t i=0;i<nst;i++)if(sc[i]>=tau){if(nnew>=ncap){ncap=ncap*5/4+1024;newS=realloc(newS,ncap*sizeof(St));}newS[nnew++]=states[i];}
        fprintf(stderr,"  kept %zu of S\n",nnew);}
      for(int chunk=0;chunk<P;chunk++){
        fclear();fcount=0;
        for(size_t i=0;i<nst;i++){int n=succ(&states[i],WMAX,out);
          for(int q=0;q<n;q++){if((int)((hsh(&out[q])>>56)%P)!=chunk)continue;if(lookup(&out[q])>=0)continue;fget(&out[q])->aest+=aT[i];}}
        for(size_t k=0;k<fcap;k++){if(ftab[k].s.r==255)continue;
          int n=succ(&ftab[k].s,WMAX,out2);double bs=0;for(int q=0;q<n;q++){int64_t j=lookup(&out2[q]);if(j>=0)bs+=bT[j];}
          float score=(float)((ftab[k].aest/sr)*(bs/sr));
          if(pass==0){nfs++;if(score>0)hist[BIN(score)]++;}
          else if(score>=tau){if(nnew>=ncap){ncap=ncap*5/4;newS=realloc(newS,ncap*sizeof(St));}newS[nnew++]=ftab[k].s;}}
        fprintf(stderr,"  pass %d chunk %d/%d: frontier entries %zu\n",pass,chunk,P,fcount);}
      if(pass==1){
        // rebuild S = newS; warm start a: old a for old states, aest-based for new (approx: 0.0 -> let iteration fill)
        // keep old a values by lookup before replacing
        float*na=malloc(nnew*4),*nb=malloc(nnew*4);
        for(size_t i=0;i<nnew;i++){int64_t j=lookup(&newS[i]);if(j>=0){na[i]=a[j];nb[i]=b[j];}else{na[i]=0;nb[i]=0;}}
        free(states);states=newS;nst=nnew;cap=ncap;
        if(nst*10>hcap*7){fprintf(stderr,"S hash too small\n");return 1;}
        hash_insert_all();build_csr();
        free(a);free(b);free(u);a=na;b=nb;u=malloc(nst*4);
        // give new states a small positive start so they are not stuck at 0 (they get filled by T^T anyway)
        for(size_t i=0;i<nst;i++)if(a[i]==0)a[i]=1e-6f;for(size_t i=0;i<nst;i++)if(b[i]==0)b[i]=1e-6f;
      }}
    free(sc);free(aT);free(bT);
    fprintf(stderr,"  stage %d -> %zu states, %zu edges [%.0fs]\n",stage,nst,nedges,now()-t0);
  }
  // write certificate: states + integer vector
  {long long pn,pd;int64_t*v;double c=certify(a,&pn,&pd,&v);
   char fn[512];snprintf(fn,512,"%s_states.bin",outp);FILE*f=fopen(fn,"wb");fwrite(states,SB,nst,f);fclose(f);
   snprintf(fn,512,"%s_vec.bin",outp);f=fopen(fn,"wb");fwrite(v,8,nst,f);fclose(f);
   printf("FINAL: states=%zu CERTIFIED rho(T^2) >= %lld/%lld = %.8f => R >= %.6f (written %s_*.bin)\n",nst,pn,pd,c,c*(1-1e-15),outp);free(v);}
  return 0;}
