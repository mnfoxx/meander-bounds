// semi-meanders (stamp 1 on top) with pruning: nonhairpins<=K or depth<=D
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
int n,K=1000,D=1000; int pile[64],pos[64]; long long cnt[64];
int crosses(int a,int b,int c,int d){return (a<c&&c<b&&b<d)||(c<a&&a<d&&d<b);}
int depthof(int k){ // max nesting depth among folds (i,i+1), i<k
  int best=0;
  for(int i=1;i<k;i++){int a=pos[i],b=pos[i+1];if(a>b){int t=a;a=b;b=t;}int c=0;
    for(int j=(i&1)?1:2;j<k;j+=2){int e=pos[j],f=pos[j+1];if(e>f){int t=e;e=f;f=t;}if(e<a&&b<f)c++;}
    if(c>best)best=c;} return best;}
void rec(int k,int nonhair){
  if(k==n){ int h=0; for(int i=1;i<n;i++) if(abs(pos[i]-pos[i+1])==1)h++; cnt[h]++; return;}
  int s=k+1;
  for(int g=1;g<=k;g++){
    for(int i=k;i>g;i--)pile[i]=pile[i-1]; pile[g]=s;
    for(int i=0;i<=k;i++)pos[pile[i]]=i;
    int a=pos[k],b=pos[s];if(a>b){int t=a;a=b;b=t;} int ok=1;
    for(int j=k-2;j>=1;j-=2){int c=pos[j],d=pos[j+1];if(c>d){int t=c;c=d;d=t;}if(crosses(a,b,c,d)){ok=0;break;}}
    if(ok){
      int nh=nonhair;
      // folds that are now permanently non-hairpin: any fold (i,i+1) with s strictly between; count only newly broken
      // simpler: recompute count of folds among 1..k-1 with |pos diff|>1 ... but that includes temporary ones. Permanent iff a stamp lies between.
      // all stamps between are permanent, so recompute: nonhair = #folds i<k with |pos(i)-pos(i+1)|>1
      nh=0; for(int i=1;i<k;i++) if(abs(pos[i]-pos[i+1])>1) nh++;
      if(nh<=K && (D>=60 || depthof(s)<=D)) rec(k+1,nh);
    }
    for(int i=g;i<k;i++)pile[i]=pile[i+1];
  }
}
int main(int argc,char**argv){n=atoi(argv[1]);if(argc>2)K=atoi(argv[2]);if(argc>3)D=atoi(argv[3]);
 pile[0]=1;pos[1]=0;rec(1,0);long long t=0;for(int i=0;i<64;i++)t+=cnt[i];
 printf("%d %lld :",n,t);for(int i=0;i<n;i++)printf(" %lld",cnt[i]);printf("\n");}
