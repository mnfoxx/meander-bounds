// Width-bounded closed-meander automaton, MIRROR QUOTIENT, packed states, exact certificate.
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
int W;
typedef struct { uint8_t r,l; uint8_t p[14]; } St; // nibble-packed labels: R labels then L labels
static inline int getl(const St*s,int i){return (s->p[i>>1]>>((i&1)*4))&15;}
static void pack(St*s,const uint8_t*lab,int n){memset(s->p,0,14);for(int i=0;i<n;i++)s->p[i>>1]|=lab[i]<<((i&1)*4);}
static uint64_t hsh(const St*a){uint64_t h=1469598103934665603ULL;const uint8_t*b=(const uint8_t*)a;for(int i=0;i<16;i++){h^=b[i];h*=1099511628211ULL;}return h;}
St*states;size_t nst=0,cap;uint32_t*htab;size_t hcap;
size_t find_or_add(const St*s){uint64_t h=hsh(s);size_t i=h&(hcap-1);
  while(htab[i]!=0xffffffffu){if(memcmp(&states[htab[i]],s,16)==0)return htab[i];i=(i+1)&(hcap-1);}
  if(nst>=cap){cap=cap*5/4;states=realloc(states,cap*sizeof(St));}
  states[nst]=*s;htab[i]=nst;return nst++;}
static void canon_lab(uint8_t*lab,int n){uint8_t map[32];memset(map,0,32);int nxt=1;for(int i=0;i<n;i++){uint8_t c=lab[i];if(!map[c])map[c]=nxt++;lab[i]=map[c];}}
static void make_state(St*t,uint8_t*R,int nr,uint8_t*L,int nl){ // canonical with mirror: min of (R,L) and (L,R)
  uint8_t a[32],b[32];memcpy(a,R,nr);memcpy(a+nr,L,nl);canon_lab(a,nr+nl);
  memcpy(b,L,nl);memcpy(b+nl,R,nr);canon_lab(b,nr+nl);
  St t1,t2;t1.r=nr;t1.l=nl;pack(&t1,a,nr+nl);t2.r=nl;t2.l=nr;pack(&t2,b,nr+nl);
  *t=(memcmp(&t1,&t2,16)<=0)?t1:t2;}
uint32_t*eoff,*edst;size_t nedges=0,ecap;
int main(int argc,char**argv){
  W=atoi(argv[1]);int iters=argc>2?atoi(argv[2]):300;
  cap=58000000;states=malloc(cap*sizeof(St));hcap=1ULL<<27;htab=malloc(hcap*4);memset(htab,0xff,hcap*4);
  ecap=165000000;edst=malloc(ecap*4);
  St init;memset(&init,0,16);find_or_add(&init);
  size_t eoffcap=58000000;eoff=malloc(eoffcap*4);
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
      if(nst*2>hcap){fprintf(stderr,"hash too small\n");return 1;}
      size_t j=find_or_add(&t);
      if(nedges>=ecap){ecap=ecap*5/4;edst=realloc(edst,ecap*4);}
      edst[nedges++]=j;}}
  eoff[nst]=nedges;free(htab);free(states);
  fprintf(stderr,"W=%d orbit-states=%zu edges=%zu\n",W,nst,nedges);
  double*v=malloc(nst*8),*u=malloc(nst*8);
  for(size_t i=0;i<nst;i++)v[i]=1.0;
  const char*ck=argc>3?argv[3]:NULL; if(ck){FILE*f=fopen(ck,"rb");if(f){size_t got=fread(v,8,nst,f);fclose(f);fprintf(stderr,"resumed %zu\n",got);}}
  double lam=1;
  for(int it=0;it<iters;it++){
    memset(u,0,nst*8);for(size_t i=0;i<nst;i++){double vi=v[i];for(uint32_t k=eoff[i];k<eoff[i+1];k++)u[edst[k]]+=vi;}
    memset(v,0,nst*8);for(size_t i=0;i<nst;i++){double ui=u[i];for(uint32_t k=eoff[i];k<eoff[i+1];k++)v[edst[k]]+=ui;}
    double mx=0;for(size_t i=0;i<nst;i++)if(v[i]>mx)mx=v[i];lam=mx;for(size_t i=0;i<nst;i++)v[i]/=mx;
    if(it%50==49)fprintf(stderr,"iter %d lam %.9f\n",it+1,lam);}
  if(ck){FILE*f=fopen(ck,"wb");fwrite(v,8,nst,f);fclose(f);fprintf(stderr,"saved\n");}
  int64_t*vi=(int64_t*)u;
  for(size_t i=0;i<nst;i++){double q=v[i]*(double)(1LL<<52);vi[i]=(int64_t)floor(q);if(vi[i]<0)vi[i]=0;}
  free(v);
  int64_t*t1=malloc(nst*8),*t2=malloc(nst*8);
  memset(t1,0,nst*8);for(size_t i=0;i<nst;i++){int64_t x=vi[i];if(x)for(uint32_t k=eoff[i];k<eoff[i+1];k++)t1[edst[k]]+=x;}
  memset(t2,0,nst*8);for(size_t i=0;i<nst;i++){int64_t x=t1[i];if(x)for(uint32_t k=eoff[i];k<eoff[i+1];k++)t2[edst[k]]+=x;}
  int64_t bn=0,bd=1;int first=1;size_t npos=0;
  for(size_t i=0;i<nst;i++){if(vi[i]<=0)continue;npos++;
    if(first){bn=t2[i];bd=vi[i];first=0;}else if((__int128)t2[i]*bd<(__int128)bn*vi[i]){bn=t2[i];bd=vi[i];}}
  double ratio=(double)bn/(double)bd;
  printf("W=%d orbit-states=%zu  rho_est=%.8f  CERTIFIED: rho(T^2) >= %lld/%lld = %.9f  => R >= %.6f  (positive: %zu of %zu)\n",
    W,nst,lam,(long long)bn,(long long)bd,ratio,ratio*(1-1e-15),npos,nst);
  return 0;}
