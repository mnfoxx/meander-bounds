#include <stdio.h>
#include <stdlib.h>
#include <string.h>
int n; int pile[64],pos[64]; long long total,bottom,bydepth[64]; double sumdepth;
int crosses(int a,int b,int c,int d){return (a<c&&c<b&&b<d)||(c<a&&a<d&&d<b);}
void stats(){ // max stack height on either side = depth+1
  for(int i=0;i<n;i++)pos[pile[i]]=i;
  int best=0;
  for(int i=1;i<n;i++){int a=pos[i],b=pos[i+1];if(a>b){int t=a;a=b;b=t;}int c=0;
    for(int j=(i&1)?1:2;j<n;j+=2){int e=pos[j],f=pos[j+1];if(e>f){int t=e;e=f;f=t;}if(e<a&&b<f)c++;}
    if(c>best)best=c;}
  total++; bydepth[best]++; sumdepth+=best; if(pos[n]==n-1) bottom++;
}
void rec(int k){
  if(k==n){stats();return;}
  int s=k+1;for(int g=1;g<=k;g++){for(int i=k;i>g;i--)pile[i]=pile[i-1];pile[g]=s;for(int i=0;i<=k;i++)pos[pile[i]]=i;
  int a=pos[k],b=pos[s];if(a>b){int t=a;a=b;b=t;}int ok=1;
  for(int j=k-2;j>=1;j-=2){int c=pos[j],d=pos[j+1];if(c>d){int t=c;c=d;d=t;}if(crosses(a,b,c,d)){ok=0;break;}}
  if(ok)rec(k+1);for(int i=g;i<k;i++)pile[i]=pile[i+1];}}
int main(int c,char**v){n=atoi(v[1]);pile[0]=1;pos[1]=0;rec(1);
 printf("%d total=%lld bottom=%lld meandepth=%.4f depthdist:",n,total,bottom,sumdepth/total);
 for(int i=0;i<n;i++)printf(" %lld",bydepth[i]);printf("\n");}
