// Width-bounded closed-meander automaton; power iteration; EXACT integer Collatz-Wielandt certificate.
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#define SMAX 30
int W;
typedef struct { uint8_t r,l; uint8_t s[SMAX]; } St;
static uint64_t hsh(const St*a){uint64_t h=1469598103934665603ULL;h^=a->r;h*=1099511628211ULL;h^=a->l;h*=1099511628211ULL;
  for(int i=0;i<a->r+a->l;i++){h^=a->s[i];h*=1099511628211ULL;}return h;}
static int eq(const St*a,const St*b){if(a->r!=b->r||a->l!=b->l)return 0;return memcmp(a->s,b->s,a->r+a->l)==0;}
St*states;size_t nst=0,cap;int64_t*htab;size_t hcap;
size_t find_or_add(const St*s){uint64_t h=hsh(s);size_t i=h&(hcap-1);
  while(htab[i]>=0){if(eq(&states[htab[i]],s))return htab[i];i=(i+1)&(hcap-1);}
  if(nst>=cap){cap*=2;states=realloc(states,cap*sizeof(St));}
  states[nst]=*s;htab[i]=nst;return nst++;}
static void canon(St*s){uint8_t map[64];memset(map,0,64);int nxt=1;int n=s->r+s->l;
  for(int i=0;i<n;i++){uint8_t c=s->s[i];if(!map[c])map[c]=nxt++;s->s[i]=map[c];}}
uint32_t*eoff,*edst;size_t nedges=0,ecap;
int main(int argc,char**argv){
  W=atoi(argv[1]);int iters=argc>2?atoi(argv[2]):300;
  cap=1<<20;states=malloc(cap*sizeof(St));hcap=1<<24;htab=malloc(hcap*sizeof(int64_t));memset(htab,0xff,hcap*sizeof(int64_t));
  ecap=1<<22;edst=malloc(ecap*sizeof(uint32_t));
  St init;memset(&init,0,sizeof init);find_or_add(&init);
  size_t eoffcap=1<<20;eoff=malloc(eoffcap*sizeof(uint32_t));
  for(size_t i=0;i<nst;i++){
    if(i+1>=eoffcap){eoffcap*=2;eoff=realloc(eoff,eoffcap*sizeof(uint32_t));}
    eoff[i]=nedges;St s=states[i];int r=s.r,l=s.l;
    for(int er=0;er<2;er++)for(int el=0;el<2;el++){
      if(!er&&r==0)continue;if(!el&&l==0)continue;
      int newlen=r+l+(er?1:-1)+(el?1:-1);if(newlen>W)continue;
      uint8_t R[SMAX],L[SMAX];int nr=r,nl=l;memcpy(R,s.s,r);memcpy(L,s.s+r,l);
      int touched[2],nt=0;if(!er)touched[nt++]=R[--nr];if(!el)touched[nt++]=L[--nl];
      if(nt==2&&touched[0]==touched[1])continue;
      int cid;
      if(nt==0){int mx=0;for(int k=0;k<nr;k++)if(R[k]>mx)mx=R[k];for(int k=0;k<nl;k++)if(L[k]>mx)mx=L[k];cid=mx+1;}
      else{cid=touched[0];if(nt==2&&touched[1]<cid)cid=touched[1];
        for(int k=0;k<nr;k++)for(int q=0;q<nt;q++)if(R[k]==touched[q])R[k]=cid;
        for(int k=0;k<nl;k++)for(int q=0;q<nt;q++)if(L[k]==touched[q])L[k]=cid;}
      if(er)R[nr++]=cid;if(el)L[nl++]=cid;
      St t;memset(&t,0,sizeof t);t.r=nr;t.l=nl;memcpy(t.s,R,nr);memcpy(t.s+nr,L,nl);canon(&t);
      if(nst*2>hcap){size_t nh=hcap*2;int64_t*nt2=malloc(nh*sizeof(int64_t));memset(nt2,0xff,nh*sizeof(int64_t));
        for(size_t k=0;k<nst;k++){uint64_t hh=hsh(&states[k]);size_t j=hh&(nh-1);while(nt2[j]>=0)j=(j+1)&(nh-1);nt2[j]=k;}
        free(htab);htab=nt2;hcap=nh;}
      size_t j=find_or_add(&t);
      if(nedges>=ecap){ecap*=2;edst=realloc(edst,ecap*sizeof(uint32_t));}
      edst[nedges++]=j;}}
  eoff[nst]=nedges;free(htab);free(states);
  fprintf(stderr,"W=%d states=%zu edges=%zu\n",W,nst,nedges);
  double*v=malloc(nst*sizeof(double)),*u=malloc(nst*sizeof(double));
  for(size_t i=0;i<nst;i++)v[i]=1.0;
  double lam=1;
  for(int it=0;it<iters;it++){
    memset(u,0,nst*sizeof(double));for(size_t i=0;i<nst;i++){double vi=v[i];for(uint32_t k=eoff[i];k<eoff[i+1];k++)u[edst[k]]+=vi;}
    memset(v,0,nst*sizeof(double));for(size_t i=0;i<nst;i++){double ui=u[i];for(uint32_t k=eoff[i];k<eoff[i+1];k++)v[edst[k]]+=ui;}
    double mx=0;for(size_t i=0;i<nst;i++)if(v[i]>mx)mx=v[i];lam=mx;for(size_t i=0;i<nst;i++)v[i]/=mx;}
  // exact certificate
  int64_t*vi=(int64_t*)u; // reuse
  for(size_t i=0;i<nst;i++){double q=v[i]*(double)(1LL<<52);vi[i]=(int64_t)floor(q);if(vi[i]<0)vi[i]=0;}
  int64_t*t1=malloc(nst*sizeof(int64_t)),*t2=malloc(nst*sizeof(int64_t));
  memset(t1,0,nst*sizeof(int64_t));for(size_t i=0;i<nst;i++){int64_t x=vi[i];if(x)for(uint32_t k=eoff[i];k<eoff[i+1];k++)t1[edst[k]]+=x;}
  memset(t2,0,nst*sizeof(int64_t));for(size_t i=0;i<nst;i++){int64_t x=t1[i];if(x)for(uint32_t k=eoff[i];k<eoff[i+1];k++)t2[edst[k]]+=x;}
  // min ratio t2[i]/vi[i] over vi>0, exact rational comparison via __int128
  int64_t bn=0,bd=1;int first=1;size_t npos=0;
  for(size_t i=0;i<nst;i++){if(vi[i]<=0)continue;npos++;
    if(first){bn=t2[i];bd=vi[i];first=0;}
    else if((__int128)t2[i]*bd<(__int128)bn*vi[i]){bn=t2[i];bd=vi[i];}}
  double ratio=(double)bn/(double)bd;
  printf("W=%d states=%zu  mu_est=%.8f   CERTIFIED: rho(M^2) >= %lld/%lld = %.9f  =>  mu >= %.8f,  R >= mu^2 >= %.6f   (states with v>0: %zu of %zu)\n",
    W,nst,sqrt(lam),(long long)bn,(long long)bd,ratio,sqrt(ratio)*(1-1e-15),ratio*(1-1e-15),npos,nst);
  return 0;}
