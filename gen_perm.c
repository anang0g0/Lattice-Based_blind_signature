#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <string.h>

#include "g_comb.c"

#define EN 1129
#define MAX 17

unsigned short prm[512][EN] = {0};
unsigned short t[512][N];

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
void random_shuffle(unsigned short *array, size_t size)
{
	for (size_t i = size; i > 1; --i)
	{
		size_t a = i - 1;
		size_t b = rand() % i;
		SWAP(int, array[a], array[b]);
	}
}

int equ(unsigned char a[EN], unsigned char b[EN])
{

	for (int i = 0; i < EN; i++)
	{
		if (a[i] != b[i])
			return 0;
	}
	return 1;
}

void beki(unsigned long long int c)
{
	int i, count = 0;
	unsigned char o[N];
	unsigned char W[N];
	for (i = 0; i < N; i++)
		W[i] = i;

	while (c > 0)
	{

		if (c % 2 == 0)
		{
			c = c >> 1;
			count++;
		}
		else
		{
			for (i = 0; i < N; i++)
				o[i] = W[t[count][i]];
			for (i = 0; i < N; i++)
				W[i] = o[i];
			c = c >> 1;
			count++;
		}
	}
}

void table(unsigned char a[EN])
{
	int i, count = 0;
	unsigned char w[N], v[N]; //, m[N], s[N];

	for (int i = 0; i < N; i++)
	{
		v[i] = i;
	}
	for (int i = 0; i < EN; i++)
	{
		t[0][i] = v[i];
		t[1][i] = v[a[i]];
		w[i] = prm[0][i];
		//s[i] = prm[1][i];
		//m[i] = prm[2][i];
		// printf("v=%d,",w[i]);
	}
	printf("\n");
	//        exit(1);

	for (count=2; count < 256;count++)
	{

		for (i = 0; i < N; i++)
		{
			v[i] = w[a[i]];
			// printf("%d,",v[i]);
		}

		printf("%d\n", count);

		for (i = 0; i < N; i++)
			a[i] = v[i];
		for (i = 0; i < N; i++)
			t[count][i] = v[i];
	}
}

void mkcycle()
{

	unsigned short p[EN], q[EN] = {0};

	// MAX個のサイクルを持つ置換を４つ出力
	int cnt = cycle(4, MAX);
	printf("%d\n", cnt);

	for (int j = 0; j < cnt; j++)
	{
		unsigned long long y = 1;
		for (int i = 0; i < MAX; i++)
		{
			printf("%d,", xx[j][i]);
			y *= xx[j][i];
		}
		printf("size[%d]=%llu\n", j, y);
	}
	// exit(1);

	for (int jj = 0; jj < cnt; jj++)
	{
		for (int i = 0; i < EN; i++)
			p[i] = i;
		random_shuffle(p, SIZE_OF_ARRAY(p));
		memset(q, 0, sizeof(q));
		int s = 0;
		for (int i = 0; i < MAX; i++)
		{
			printf("(");
			int k = xx[jj][i];
			for (int j = s; j < k + s; j++)
			{
				printf("%d,", p[j]);
				if (j != s + k - 1)
					q[p[j]] = p[(j + 1)];

				if (j == s + k - 1)
					q[p[j]] = p[s];
			}

			printf(")");
			s += k;
		}
		printf("\n");
		for (int i = 0; i < EN; i++)
		{
			prm[jj][i] = q[i];
		}
		for (int i = 0; i < EN; i++)
		{
			printf("%d,", prm[jj][i]);
		}
		printf("\n");
	}

	return;
}

void nike(int q)
{

	for (int i = 0; i < EN; i++)
		printf(".%d,", prm[q][prm[q][i]]);
	printf("\n");
}

int main()
{
	unsigned char o[EN];

	srand(clock());
    printf("要素の数を入力してください: ");
    scanf("%d", &n);
    permutation(0);
	l=1;
	for(int i=1;i<=n;i++)
	l*=i;
	for(int i=0;i<l;i++){
	for(int j=0;j<n;j++)
	printf("%d ",vv[i].x[j]);
	printf("\n");
	}
	//exit(1);
	int cnt=0;
	vec v={0};
	for(int i=0;i<n;i++)
	v.x[i]=i+1;
  do {
	vx[cnt++]=v;
  } while (next_permutation(0, n, v.x));
	for(int i=0;i<l;i++){
	for(int j=0;j<n;j++)
	printf("%d ",vx[i].x[j]);
	printf("\n");
	}
	exit(1);

	mkcycle();
	//nike(1);
	//for (int i = 0; i < EN; i++)
	//	o[i] = prm[1][i];
	//table(o);

	return 0;
}

