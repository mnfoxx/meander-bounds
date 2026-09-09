#include <stdio.h>
#include <stdlib.h>
int n,D; int pile[64],pos[64]; long long cnt;
int crosses(int a,int b,int c,int d){return (a<c&&c<b&&b<d)||(c<a&&a<d&&d<b);}
int depthof(int k){int best=0;for(int i=1;i<k;i++){int a=pos[i],b=pos[i+1];if(a>b){int t=a;a=b;b=t;}int c=0;
 for(int j=(i&1)?1:2;j<k;j+=2){int e=pos[j],f=pos[j+1];if(e>f){int t=e;e=f;f=t;}if(e<a&&b<f)c++;}if(c>best)best=c;}return best;}
void rec(int k){
 if(k==n){if(pos[n]==n-1)cnt++;return;}
 int s=k+1;for(int g=1;g<=k;g++){for(int i=k;i>g;i--)pile[i]=pile[i-1];pile[g]=s;for(int i=0;i<=k;i++)pos[pile[i]]=i;
 int a=pos[k],b=pos[s];if(a>b){int t=a;a=b;b=t;}int ok=1;
 for(int j=k-2;j>=1;j-=2){int c=pos[j],d=pos[j+1];if(c>d){int t=c;c=d;d=t;}if(crosses(a,b,c,d)){ok=0;break;}}
 if(ok&&depthof(s)<=D)rec(k+1);for(int i=g;i<k;i++)pile[i]=pile[i+1];}}
int main(int c,char**v){n=atoi(v[1]);D=atoi(v[2]);pile[0]=1;pos[1]=0;rec(1);printf("%lld ",cnt);}
