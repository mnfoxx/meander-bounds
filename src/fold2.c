#include <stdio.h>
#include <stdlib.h>
#include <string.h>
int n; int pile[64]; int pos[64];
long long by2[64], byhair[64], bydepth[64], bylen[4096], bybottom[64];
int crosses(int a,int b,int c,int d){return (a<c&&c<b&&b<d)||(c<a&&a<d&&d<b);}
void stats(){
  for(int i=0;i<n;i++) pos[pile[i]]=i;
  if(pos[1]!=0) return; // semi-meander: stamp 1 on top
  int hair=0,len=0,depth=0;
  for(int i=1;i<n;i++){int d=abs(pos[i]-pos[i+1]); if(d==1)hair++; len+=d;}
  // nesting depth: for each fold count folds on same side containing it
  for(int i=1;i<n;i++){int a=pos[i],b=pos[i+1];if(a>b){int t=a;a=b;b=t;} int c=0;
    for(int j=i%2?1:2;j<n;j+=2){int e=pos[j],f=pos[j+1];if(e>f){int t=e;e=f;f=t;} if(e<a&&b<f)c++;}
    if(c>depth)depth=c;}
  by2[pos[2]]++; byhair[hair]++; bydepth[depth]++; bylen[len]++; bybottom[pile[n-1]]++;
}
void rec(int k){
  if(k==n){stats();return;}
  int s=k+1;
  for(int g=0;g<=k;g++){
    if(g==0 && k>=1) { /* stamp 1 must stay on top: only allow g>=1 */ continue; }
    for(int i=k;i>g;i--)pile[i]=pile[i-1]; pile[g]=s;
    for(int i=0;i<=k;i++)pos[pile[i]]=i;
    int a=pos[k],b=pos[s];if(a>b){int t=a;a=b;b=t;} int ok=1;
    for(int j=k-2;j>=1;j-=2){int c=pos[j],d=pos[j+1];if(c>d){int t=c;c=d;d=t;} if(crosses(a,b,c,d)){ok=0;break;}}
    if(ok)rec(k+1);
    for(int i=g;i<k;i++)pile[i]=pile[i+1];
  }
}
int main(int argc,char**argv){n=atoi(argv[1]);pile[0]=1;rec(1);
 long long t=0; for(int i=0;i<n;i++)t+=by2[i];
 printf("n=%d semimeanders=%lld\n",n,t);
 printf("pos2:"); for(int i=1;i<n;i++)printf(" %lld",by2[i]); printf("\n");
 printf("bottom:"); for(int i=1;i<=n;i++)printf(" %lld",bybottom[i]); printf("\n");
 printf("hairpins:"); for(int i=0;i<n;i++)printf(" %lld",byhair[i]); printf("\n");
 printf("depth:"); for(int i=0;i<n;i++)printf(" %lld",bydepth[i]); printf("\n");
 printf("len:"); for(int i=0;i<4096;i++)if(bylen[i])printf(" %d:%lld",i,bylen[i]); printf("\n\n");
}
