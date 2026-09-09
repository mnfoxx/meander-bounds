// Rootless bounded-height meander automaton in C. States: two stacks (R,L) of path labels, canonical.
// Build by BFS with hash table; then power iteration on M^2 (bipartite); report Collatz-Wielandt lower bound.
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
int H; // max stack height
typedef struct { uint8_t r,l; uint8_t s[40]; } St; // r = R height, l = L height, s = R labels then L labels
static uint64_t hsh(const St*a){uint64_t h=1469598103934665603ULL; h^=a->r; h*=1099511628211ULL; h^=a->l; h*=1099511628211ULL;
  for(int i=0;i<a->r+a->l;i++){h^=a->s[i];h*=1099511628211ULL;} return h;}
static int eq(const St*a,const St*b){ if(a->r!=b->r||a->l!=b->l) return 0; return memcmp(a->s,b->s,a->r+a->l)==0; }
St *states; size_t nst=0, cap;
int64_t *htab; size_t hcap;
size_t find_or_add(const St*s){ uint64_t h=hsh(s); size_t i=h&(hcap-1);
  while(htab[i]>=0){ if(eq(&states[htab[i]],s)) return htab[i]; i=(i+1)&(hcap-1);} 
  if(nst>=cap){cap*=2; states=realloc(states,cap*sizeof(St));}
  states[nst]=*s; htab[i]=nst; return nst++; }
static void canon(St*s){ // relabel by first appearance
  uint8_t map[64]; memset(map,0,64); int nxt=1; int n=s->r+s->l;
  for(int i=0;i<n;i++){ uint8_t c=s->s[i]; if(!map[c]) map[c]=nxt++; s->s[i]=map[c]; } }
// edges: CSR
uint32_t *eoff,*edst; size_t nedges=0, ecap;
int main(int argc,char**argv){
  H=atoi(argv[1]); int iters=argc>2?atoi(argv[2]):200;
  cap=1<<20; states=malloc(cap*sizeof(St)); hcap=1<<24; htab=malloc(hcap*sizeof(int64_t)); memset(htab,0xff,hcap*sizeof(int64_t));
  ecap=1<<22; edst=malloc(ecap*sizeof(uint32_t)); 
  St init; memset(&init,0,sizeof init); find_or_add(&init);
  size_t eoffcap=1<<20; eoff=malloc(eoffcap*sizeof(uint32_t));
  for(size_t i=0;i<nst;i++){
    if(i+1>=eoffcap){eoffcap*=2;eoff=realloc(eoff,eoffcap*sizeof(uint32_t));}
    eoff[i]=nedges;
    St s=states[i]; int r=s.r,l=s.l;
    for(int er=0;er<2;er++)for(int el=0;el<2;el++){ // 1=open,0=close
      if(!er && r==0) continue; if(!el && l==0) continue;
      if(er && r+l+(el?1:0)>H) continue; if(el && r+l+(er?1:0)>H) continue; if(er&&el&&r+l+2>H) continue;
      St t; memset(&t,0,sizeof t);
      uint8_t R[40],L[40]; int nr=r,nl=l; memcpy(R,s.s,r); memcpy(L,s.s+r,l);
      int touched[2],nt=0;
      if(!er) touched[nt++]=R[--nr];
      if(!el) touched[nt++]=L[--nl];
      if(nt==2 && touched[0]==touched[1]) continue; // cycle
      int cid;
      if(nt==0){ int mx=0; for(int k=0;k<nr;k++) if(R[k]>mx)mx=R[k]; for(int k=0;k<nl;k++) if(L[k]>mx)mx=L[k]; cid=mx+1; }
      else { cid=touched[0]; if(nt==2 && touched[1]<cid) cid=touched[1];
        for(int k=0;k<nr;k++) for(int q=0;q<nt;q++) if(R[k]==touched[q]) R[k]=cid;
        for(int k=0;k<nl;k++) for(int q=0;q<nt;q++) if(L[k]==touched[q]) L[k]=cid; }
      if(er) R[nr++]=cid; if(el) L[nl++]=cid;
      t.r=nr;t.l=nl; memcpy(t.s,R,nr); memcpy(t.s+nr,L,nl); canon(&t);
      // rehash table grow if needed
      if(nst*2>hcap){ size_t nh=hcap*2; int64_t*nt2=malloc(nh*sizeof(int64_t)); memset(nt2,0xff,nh*sizeof(int64_t));
        for(size_t k=0;k<nst;k++){uint64_t hh=hsh(&states[k]); size_t j=hh&(nh-1); while(nt2[j]>=0) j=(j+1)&(nh-1); nt2[j]=k;}
        free(htab); htab=nt2; hcap=nh; }
      size_t j=find_or_add(&t);
      if(nedges>=ecap){ecap*=2;edst=realloc(edst,ecap*sizeof(uint32_t));}
      edst[nedges++]=j;
    }
  }
  eoff[nst]=nedges;
  fprintf(stderr,"H=%d states=%zu edges=%zu\n",H,nst,nedges);
  // power iteration on M^2 (exclude the transient empty state: it has no incoming edges, fine)
  double *v=malloc(nst*sizeof(double)),*w=malloc(nst*sizeof(double)),*u=malloc(nst*sizeof(double));
  for(size_t i=0;i<nst;i++)v[i]=1.0;
  double lam=0;
  for(int it=0;it<iters;it++){
    memset(u,0,nst*sizeof(double)); for(size_t i=0;i<nst;i++) for(uint32_t k=eoff[i];k<eoff[i+1];k++) u[edst[k]]+=v[i];
    memset(w,0,nst*sizeof(double)); for(size_t i=0;i<nst;i++) for(uint32_t k=eoff[i];k<eoff[i+1];k++) w[edst[k]]+=u[i];
    double nrm=0; for(size_t i=0;i<nst;i++) nrm+=w[i]*w[i]; nrm=sqrt(nrm); lam=nrm; for(size_t i=0;i<nst;i++) v[i]=w[i]/nrm;
  }
  // Collatz-Wielandt lower bound on rho(M^2) using v (v>0 on recurrent states): min over i with v_i>0 of (M^2 v)_i / v_i
  memset(u,0,nst*sizeof(double)); for(size_t i=0;i<nst;i++) for(uint32_t k=eoff[i];k<eoff[i+1];k++) u[edst[k]]+=v[i];
  memset(w,0,nst*sizeof(double)); for(size_t i=0;i<nst;i++) for(uint32_t k=eoff[i];k<eoff[i+1];k++) w[edst[k]]+=u[i];
  // note: above computes M^T-style propagation (v_i pushed to successors) => w = (M^T)^2 v; CW bound works for M^T too (same spectrum)
  double cw=1e300, cwmax=0; size_t npos=0;
  for(size_t i=1;i<nst;i++){ if(v[i]>1e-12){ double q=w[i]/v[i]; if(q<cw)cw=q; if(q>cwmax)cwmax=q; npos++; } }
  printf("D=%d states=%zu  mu_est=%.8f  CW lower bound on mu: %.8f  (upper %.8f)  positive states=%zu\n",H-1,nst,sqrt(lam),sqrt(cw),sqrt(cwmax),npos);
  return 0;
}
