#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "sha3.c"

#define M 128
#define K 64
#define N 106
#define E 7

// #include "aaa.c"

#define MATRIX_SIZE 128*64
#define SHM_KEY 1234

typedef struct
{
    short x[N][N];
} MTX2;

#define SIZE_OF_ARRAY(array) (sizeof(array) / sizeof(array[0]))
#define SWAP(type, a, b) \
    {                    \
        type work = a;   \
        a = b;           \
        b = work;        \
    }

/*
    Fisher-Yates shuffle による方法
    配列の要素をランダムシャッフルする
*/
void random_shuffle(unsigned short *array, size_t size)
{
    for (size_t i = size; i > 1; --i)
    {
        size_t a = i - 1;
        size_t b = rand() % i;
        SWAP(int, array[a], array[b]);
    }
}

long long sum(short a[N]){
    long long u=0;
    for(int i=0;i<N;i++)
        u+=(a[i]+a[(i+1)%N]);

    return u;
}

void idm(){
    short z[N]={0};
    short v[N]={0};
    short ra[N]={0};
    short w[N]={0};

    int i, j, k;
    short P[N] = {0},Pa[N]={0};
    short b = 0;
    int a2[32][64] = {0};
    int a3[32][64] = {0};


    for (i = 0; i < N; i++){
        P[i] = i;
        Pa[i]=i;
    }
    for (i = 0; i < N; i++)
        P[i] = i;    
    random_shuffle(P, N);
    random_shuffle(Pa,N);
    short inv_P[N];
    short inv_Pa[N]={0};
    for (i = 0; i < N; i++)
    {
        printf("%d,", P[i]);
        inv_P[P[i]] = i;
    }

    for(i=0;i<N;i++){
    ra[i]=rand()%4096;
    w[i]=rand()%4096;
    if(rand()%2==1)
    ra[i]= -1*ra[i];
    if(rand()%2==1)
    v[i]= -1*w[inv_Pa[i]];
    }
    for(i=0;i<N;i++)
    inv_Pa[Pa[i]]=i;
    printf("\n");


    short p1[N]={0},p2[N]={0},p3[N]={0};
    int ca;
    
    for(j=0;j<1;j++){    
    ca=rand()%4093;
    for(int kk=0;kk<1;kk++)
    {
    for(i=0;i<N;i++){
    z[i]=ra[P[i]]+ca*v[Pa[P[i]]];
    //printf("aa%d %d\n",z[i],p3[i]);
    }

    for(i=0;i<N;i++){
    p1[i]=ca*v[Pa[i]];
    p2[i]=ra[P[i]];
    p3[i]=ca*v[Pa[P[i]]];
    //z[i]=ca*v[Pa[i]];
    }
    //for(i=0;i<N;i++)
    //printf("%d %d\n",p3[i],p2[i]);
    long long c0=sum(p3);
    long long c1=sum(p2);
    printf("%lld\n",c0);
    printf("%lld\n",c1);
b=rand()%2;
short tmp[N]={0};
if(b==0)
{
    for(i=0;i<N;i++)
    tmp[i]=z[inv_P[i]];
    if(sum(ra)==(sum(tmp)-sum(p1))){
    printf("That's True\n");
    }else{
        printf("baka1\n");
        exit(1);
    }
}
if(b==1){
    if(sum(p2)==(sum(z)-sum(p3)))
    {
        printf("That's True!\n");
    }else{
        printf("baka2\n");
        exit(1);
    }
}
    }
    printf("::\n");
    }


}


int main()
{
    short z[17][17][N]={0};
    short v[N]={0};
    short ra[N]={0};
    short w[N]={0};

    int i, j, k;
    short P[N] = {0},Pa[N]={0};
    short b[N][N] = {0};
    int a2[32][64] = {0};
    int a3[32][64] = {0};

    srand(clock());

    idm();
    exit(1);

    for (i = 0; i < N; i++){
        P[i] = i;
        Pa[i]=i;
    }
    for (i = 0; i < N; i++)
        P[i] = i;    
    random_shuffle(P, N);
    random_shuffle(Pa,N);
    short inv_P[N];
    short inv_Pa[N]={0};
    for (i = 0; i < N; i++)
    {
        printf("%d,", P[i]);
        inv_P[P[i]] = i;
    }

    for(i=0;i<N;i++){
    ra[i]=rand()%4096;
    w[i]=rand()%4096;
    if(rand()%2==1)
    ra[i]= -1*ra[i];
    if(rand()%2==1)
    v[i]= -1*w[inv_Pa[i]];
    }
    for(i=0;i<N;i++)
    inv_Pa[Pa[i]]=i;
    printf("\n");


    short p1[N]={0},p2[N]={0},p3[N]={0};
    int ca;
    
    for(j=0;j<17;j++){    
    ca=rand()%4093;
    for(int kk=0;kk<17;kk++)
    {
    for(i=0;i<N;i++){
    z[j][kk][i]=ra[P[i]]+ca*v[Pa[P[i]]];
    //printf("aa%d %d\n",z[i],p3[i]);
    }

    for(i=0;i<N;i++){
    p1[i]=v[Pa[i]];
    p2[i]=ra[i]+p1[i]*ca;
    p3[i]=ra[P[i]];
    //z[i]=ca*v[Pa[i]];
    }
    for(i=0;i<N;i++)
    printf("%d %d\n",p3[i],p2[i]);
    long long c0=sum(p3);
    long long c1=sum(p2);
    printf("%lld\n",c0);
    printf("%lld\n",c1);
    
    short va[N]={0};
    for(i=0;i<N;i++)
    va[i]=z[j][kk][inv_P[i]];
    printf("%lld\n",sum(va));
    short vb[N]={0};
    for(i=0;i<N;i++)
    vb[i]=z[j][kk][i]-ca*v[Pa[P[i]]];
    printf("%lld\n",sum(vb));
    a2[j][kk]=sum(va);
    a3[j][kk]=sum(vb);

    for(k=0;k<N;k++){
    ra[k]+=ra[P[k]];
    v[i]+=v[Pa[P[i]]];
    }
    }
    printf("::\n");
    }
    printf("\n");
    for(j=0;j<17;j++){
    for(i=0;i<17;i++)
    printf("%d,%d\n",a2[j][i],a3[j][i]);
    //for(int kk=0;kk<17;kk++)
    //for(i=0;i<106;i++)
    //printf("%d,",z[j][kk][i]);
    //printf("\n");
    }

    return 0;
}

