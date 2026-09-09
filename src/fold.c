// Enumerate stamp foldings by inserting stamps 1..n into a pile.
// pile[] = labels top-to-bottom. Fold between i and i+1 lies on side (i&1).
// Two same-side folds must not cross (interleave).
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
int n; int pile[64]; int pos[64]; 
long long total; long long byfirst[64]; long long bylast[64]; long long bytop[64];
int crosses(int a,int b,int c,int d){ // intervals [a,b],[c,d] with a<b,c<d
  return (a<c && c<b && b<d) || (c<a && a<d && d<b);
}
void rec(int k){ // stamps 1..k placed; k stamps in pile
  if(k==n){
    total++;
    for(int i=0;i<n;i++) pos[pile[i]]=i;
    byfirst[pos[1]]++; bylast[pos[n]]++; bytop[pile[0]]++;
    return;
  }
  int s=k+1;
  for(int g=0; g<=k; g++){
    // insert s at gap g (g = index in new pile)
    for(int i=k;i>g;i--) pile[i]=pile[i-1];
    pile[g]=s;
    for(int i=0;i<=k;i++) pos[pile[i]]=i;
    // new fold (k, s), side k&1
    int a=pos[k], b=pos[s]; if(a>b){int t=a;a=b;b=t;}
    int ok=1;
    for(int j=k-2; j>=1; j-=2){ // same side folds (j,j+1)
      int c=pos[j], d=pos[j+1]; if(c>d){int t=c;c=d;d=t;}
      if(crosses(a,b,c,d)){ok=0;break;}
    }
    if(ok) rec(k+1);
    for(int i=g;i<k;i++) pile[i]=pile[i+1];
  }
}
int main(int argc,char**argv){
  n=atoi(argv[1]); total=0; memset(byfirst,0,sizeof byfirst);
  pile[0]=1; rec(1);
  printf("n=%d total=%lld\n",n,total);
  printf("pos_of_stamp1:"); for(int i=0;i<n;i++) printf(" %lld",byfirst[i]); printf("\n");
  printf("pos_of_stampn:"); for(int i=0;i<n;i++) printf(" %lld",bylast[i]); printf("\n");
  printf("top_stamp:    "); for(int i=1;i<=n;i++) printf(" %lld",bytop[i]); printf("\n");
}
