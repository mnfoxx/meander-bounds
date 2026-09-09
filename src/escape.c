// For each semi-meander of n stamps (stamp 1 on top), compute the escape number:
// minimal k such that appending stamps n+1..n+k yields stamp n+k at the bottom (1 stays on top).
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define MAXL 40
int n; int pile[MAXL],pos[MAXL];
long long total, hist[MAXL]; double sumesc;
int crosses(int a,int b,int c,int d){return (a<c&&c<b&&b<d)||(c<a&&a<d&&d<b);}
// legal insertion of stamp L+1 at gap g (0..len) into pile p of length len; returns 1 if legal
int legal(int *p,int len,int g){
  int q[MAXL],ps[MAXL+1]; int L=len; // new stamp label L+1
  for(int i=0;i<g;i++)q[i]=p[i]; q[g]=L+1; for(int i=g;i<len;i++)q[i+1]=p[i];
  for(int i=0;i<=len;i++)ps[q[i]]=i;
  int a=ps[L],b=ps[L+1];if(a>b){int t=a;a=b;b=t;}
  for(int j=L-2;j>=1;j-=2){int c=ps[j],d=ps[j+1];if(c>d){int t=c;c=d;d=t;}if(crosses(a,b,c,d))return 0;}
  return 1;
}
// BFS over piles
typedef struct{int len; int p[MAXL];} St;
St *cur,*nxt; int ncur,nnxt; 
// visited: simple hash set of strings
#define HS (1<<22)
unsigned long long *hkeys; 
unsigned long long hsh(St*s){unsigned long long h=1469598103934665603ULL;for(int i=0;i<s->len;i++){h^=s->p[i];h*=1099511628211ULL;}return h|1;}
int seen(St*s){unsigned long long h=hsh(s);unsigned idx=h&(HS-1);while(hkeys[idx]){if(hkeys[idx]==h)return 1;idx=(idx+1)&(HS-1);}hkeys[idx]=h;return 0;}
int escape(){
  ncur=1; cur[0].len=n; memcpy(cur[0].p,pile,n*sizeof(int)); memset(hkeys,0,HS*sizeof(unsigned long long)); seen(&cur[0]);
  for(int k=0;;k++){
    nnxt=0;
    for(int i=0;i<ncur;i++){St*s=&cur[i]; if(s->p[s->len-1]==s->len) return k; // last stamp at bottom
      if(s->len+1>=MAXL) continue;
      for(int g=1;g<=s->len;g++) if(legal(s->p,s->len,g)){
        St t; t.len=s->len+1; for(int j=0;j<g;j++)t.p[j]=s->p[j]; t.p[g]=s->len+1; for(int j=g;j<s->len;j++)t.p[j+1]=s->p[j];
        if(!seen(&t)){ if(nnxt<2000000) nxt[nnxt++]=t; }
      }}
    if(nnxt==0) return -1;
    St*tmp=cur;cur=nxt;nxt=tmp;ncur=nnxt;
  }
}
void rec(int k){
  if(k==n){ for(int i=0;i<n;i++)pos[pile[i]]=i; int e=escape(); total++; hist[e]++; sumesc+=e; return;}
  int s=k+1;for(int g=1;g<=k;g++){for(int i=k;i>g;i--)pile[i]=pile[i-1];pile[g]=s;for(int i=0;i<=k;i++)pos[pile[i]]=i;
  int a=pos[k],b=pos[s];if(a>b){int t=a;a=b;b=t;}int ok=1;
  for(int j=k-2;j>=1;j-=2){int c=pos[j],d=pos[j+1];if(c>d){int t=c;c=d;d=t;}if(crosses(a,b,c,d)){ok=0;break;}}
  if(ok)rec(k+1);for(int i=g;i<k;i++)pile[i]=pile[i+1];}}
int main(int c,char**v){n=atoi(v[1]);cur=malloc(2000000*sizeof(St));nxt=malloc(2000000*sizeof(St));hkeys=malloc(HS*sizeof(unsigned long long));
 pile[0]=1;pos[1]=0;rec(1);
 printf("%d total=%lld mean_escape=%.4f hist:",n,total,sumesc/total);for(int i=0;i<n;i++)printf(" %lld",hist[i]);printf("\n");}
