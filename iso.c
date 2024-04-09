#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

#define N 5

#define SIZE_OF_ARRAY(array) (sizeof(array) / sizeof(short))
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
void random_shuffle(unsigned char *array, size_t size)
{
    for (size_t i = size; i > 1; --i)
    {
        size_t a = i - 1;
        size_t b = rand() % i;
        SWAP(int, array[a], array[b]);
    }
}

int chk(unsigned char  *a){
    for(int i=0;i<N;i++){
    if(a[i]==i)
    return -1;
    }

    return 0;
}

typedef struct {
    unsigned char v[N];
} vec;
 
typedef struct {
    unsigned char x[N];
    unsigned char pi[N];
} comb;

unsigned char rotl(unsigned char x, unsigned char r)
{
	if (r == 0)
		return x;
	if (r < 0)
		return (x >> r) | (x << (8 - r));

	return (x << r) | (x >> (8 - r));
}

comb xxx(comb a, comb b){
comb t;
int i,j;
comb inv_b={0},s={0};

for(i=0;i<N;i++)
inv_b.pi[b.pi[i]]=i;
for(i=0;i<N;i++)
s.pi[i]=b.pi[a.pi[inv_b.pi[i]]];
memcpy(a.pi,s.pi,N);
t=a;

for(i=0;i<N;i++)
t.x[i]^=rotl(b.x[i]+a.x[a.pi[i]],i%8);

return t;
}

FILE *fp;
void printx(comb a){
    //for(int i=0;i<N;i++)
    //printf("%d,",a.x[i]);
    //printf("\n");
    fwrite(a.x,1,N,fp);
    //for(int i=0;i<N;i++)
    //printf("%d,",a.pi[i]);
    //printf("\n");
}

int main(){
    int i,j;
    comb x,p,q,inv_p;
    comb t;


    fp=fopen("test.bin","wb");
    srand(clock());
    for(i=0;i<N;i++){
    p.x[i]=rand()%2;
    p.pi[i]=i;
    q.x[i]=rand()%2;
    q.pi[i]=i;
    }
    j=-1;
    while(j<0){
    random_shuffle(p.pi,N);
    j=chk(p.pi);
    printx(p);
    }
    j=-1;
    while(j<0){
    random_shuffle(q.pi,N);
    j=chk(q.pi);
    printx(q);
    }

    i=0;
    while(i<2000000){
    p=xxx(p,q);
    printx(p);
    //q=xxx(q,p);
    //printx(q);
    i++;
    }
    fclose(fp);

return 0;
}