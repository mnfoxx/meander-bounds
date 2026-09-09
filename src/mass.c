#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#define SMAX 30
int W;
typedef struct{uint8_t r,l;uint8_t s[SMAX];}St;
static uint64_t hsh(const St*a){uint64_t h=1469598103934665603ULL;h^=a->r;h*=1099511628211ULL;h^=a->l;h*=1099511628211ULL;for(int i=0;i<a->r+a->l;i++){h^=a->s[i];h*=1099511628211ULL;}return h;}
static int eq(const St*a,const St*b){if(a->r!=b->r||a->l!=b->l)return 0;return memcmp(a->s,b->s,a->r+a->l)==0;}
St*states;size_t nst=0,cap;int64_t*htab;size_t hcap;
size_t find_or_add(const St*s){uint64_t h=hsh(s);size_t i=h&(hcap-1);while(htab[i]>=0){if(eq(&states[htab[i]],s))return htab[i];i=(i+1)&(hcap-1);}if(nst>=cap){cap*=2;states=realloc(states,cap*sizeof(St));}states[nst]=*s;htab[i]=nst;return nst++;}
static void canon(St*s){uint8_t map[64];memset(map,0,64);int nxt=1;int n=s->r+s->l;for(int i=0;i<n;i++){uint8_t c=s->s[i];if(!map[c])map[c]=nxt++;s->s[i]=map[c];}}
uint32_t*eoff,*edst;size_t nedges=0,ecap;
int main(int argc,char**argv){W=atoi(argv[1]);int iters=atoi(argv[2]);
 cap=1<<20;states=malloc(cap*sizeof(St));hcap=1<<25;htab=malloc(hcap*8);memset(htab,0xff,hcap*8);ecap=1<<22;edst=malloc(ecap*4);
 St init;memset(&init,0,sizeof init);find_or_add(&init);size_t eoffcap=1<<20;eoff=malloc(eoffcap*4);
 for(size_t i=0;i<nst;i++){if(i+1>=eoffcap){eoffcap*=2;eoff=realloc(eoff,eoffcap*4);}eoff[i]=nedges;St s=states[i];int r=s.r,l=s.l;
  for(int er=0;er<2;er++)for(int el=0;el<2;el++){if(!er&&r==0)continue;if(!el&&l==0)continue;if(r+l+(er?1:-1)+(el?1:-1)>W)continue;
   uint8_t R[SMAX],L[SMAX];int nr=r,nl=l;memcpy(R,s.s,r);memcpy(L,s.s+r,l);int touched[2],nt=0;if(!er)touched[nt++]=R[--nr];if(!el)touched[nt++]=L[--nl];
   if(nt==2&&touched[0]==touched[1])continue;int cid;
   if(nt==0){int mx=0;for(int k=0;k<nr;k++)if(R[k]>mx)mx=R[k];for(int k=0;k<nl;k++)if(L[k]>mx)mx=L[k];cid=mx+1;}
   else{cid=touched[0];if(nt==2&&touched[1]<cid)cid=touched[1];for(int k=0;k<nr;k++)for(int q=0;q<nt;q++)if(R[k]==touched[q])R[k]=cid;for(int k=0;k<nl;k++)for(int q=0;q<nt;q++)if(L[k]==touched[q])L[k]=cid;}
   if(er)R[nr++]=cid;if(el)L[nl++]=cid;St t;memset(&t,0,sizeof t);t.r=nr;t.l=nl;memcpy(t.s,R,nr);memcpy(t.s+nr,L,nl);canon(&t);
   size_t j=find_or_add(&t);if(nedges>=ecap){ecap*=2;edst=realloc(edst,ecap*4);}edst[nedges++]=j;}}
 eoff[nst]=nedges;
 double*a=malloc(nst*8),*b=malloc(nst*8),*u=malloc(nst*8);for(size_t i=0;i<nst;i++){a[i]=1;b[i]=1;}
 for(int it=0;it<iters;it++){ // a: left (push forward), b: right (pull back)
  memset(u,0,nst*8);for(size_t i=0;i<nst;i++)for(uint32_t k=eoff[i];k<eoff[i+1];k++)u[edst[k]]+=a[i];double m=0;for(size_t i=0;i<nst;i++)if(u[i]>m)m=u[i];for(size_t i=0;i<nst;i++)a[i]=u[i]/m;
  memset(u,0,nst*8);for(size_t i=0;i<nst;i++){double sum=0;for(uint32_t k=eoff[i];k<eoff[i+1];k++)sum+=b[edst[k]];u[i]=sum;}m=0;for(size_t i=0;i<nst;i++)if(u[i]>m)m=u[i];for(size_t i=0;i<nst;i++)b[i]=u[i]/m;}
 double tot=0,bypaths[40]={0},bywidth[40]={0},cntpaths[40]={0};
 for(size_t i=0;i<nst;i++){St*s=&states[i];int mx=0;for(int k=0;k<s->r+s->l;k++)if(s->s[k]>mx)mx=s->s[k];double m=a[i]*b[i];tot+=m;bypaths[mx]+=m;bywidth[s->r+s->l]+=m;cntpaths[mx]+=1;}
 printf("W=%d states=%zu\n by #paths (mass%%, states):",W,nst);for(int p=0;p<=W/2;p++)printf(" %d:%.2f%%/%.0f",p,100*bypaths[p]/tot,cntpaths[p]);printf("\n by width:");for(int w=0;w<=W;w++)if(bywidth[w]>0)printf(" %d:%.1f%%",w,100*bywidth[w]/tot);printf("\n");}
