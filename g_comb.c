/*
 * Author: Hiroyuki Chishiro
 * License: 2-Clause BSD
 */
#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>

// modified. use ./a.out|grep @
#define N 1129

unsigned short xx[256][N] = {0};

typedef struct {
	unsigned short x[N];
} vec;

void swap(unsigned short *pa, unsigned short *pb)
{
	int tmp;

	tmp = *pa;
	*pa = *pb;
	*pb = tmp;
}

void rotate(size_t first, size_t middle, size_t last, unsigned short v[])
{
	size_t middle_org;

	if (first == middle || middle == last)
	{
		return;
	}

	middle_org = middle;

	while (first != middle_org && middle != last)
	{
		swap(&v[first++], &v[middle++]);
	}

	if (first == middle_org)
	{
		rotate(first, middle, last, v);
	}
	else
	{
		rotate(first, middle_org, last, v);
	}
}

bool next_combination(size_t first, size_t last, size_t r, unsigned short v[])
{
	size_t subset = first + r;
	size_t src = subset;
	size_t dst = subset;

	if (first == last || first == subset || last == subset)
	{
		return false;
	}

	while (first != src)
	{
		src--;

		if (v[src] < v[last - 1])
		{

			while (v[src] >= v[dst])
			{
				dst++;
			}

			swap(&v[src], &v[dst]);
			rotate(src + 1, dst + 1, last, v);
			rotate(subset, subset + (last - dst) - 1, last, v);

			return true;
		}
	}

	rotate(first, subset, last, v);

	return false;
}

int cycle(unsigned int rr, size_t r)
{
	unsigned short v[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193};

	size_t i;
	size_t n = 29; // sizeof(v) / sizeof(v[0]);
	// size_t r=; // = 3;
	int j = 0, count = 0;
	// for(r=2;r<16;r++)

	count = 0;
	do
	{
		j = 0;
		for (i = 0; i < r; i++)
		{
			j += v[i];
		}
		if (j == N)
		{
			for (i = 0; i < r; i++)
			{
				printf("%d ", v[i]);
				xx[count][i] = v[i];
			}
			count++;
			printf("@=%d rr=%d\n", N, rr);
			if (count == rr)
				break;
		}
		// printf("\n");
	} while (next_combination(0, n, r, v));
	printf("@count=%d r=%lu\n", count, r);
	if (rr != count)
		exit(1);

	return count;
}

unsigned long long l=0;
#define MAX_N N
int n; // 配列の要素数
vec arr; // 配列
int used[MAX_N]; // 選ばれた要素をマークする配列
int count=0;
vec vv[362880]={0};
vec vx[362880]={0};
void permutation(int depth) {
    if (depth == n) {
		vv[count]=arr;
		count++;
        return;
    }
    for (int i = 0; i < n; i++) {
        if (!used[i]) {
            arr.x[depth] = i + 1;
            used[i] = 1;
            permutation(depth + 1);
            used[i] = 0;
        }
    }
}


void reverse(size_t first, size_t last, unsigned short v[])
{
  while (first != last && first != --last) {
    swap(&v[first], &v[last]);
    first++;
  }
}

bool next_permutation(size_t first, size_t last, unsigned short v[])
{
  size_t i, j, k;

  if (first == last) {
    return false;
  }

  if (first + 1 == last) {
    return false;
  }

  i = last - 1;

  while (true) {
    j = i--;

    if (v[i] < v[j]) {
      k = last;

      while (!(v[i] < v[--k])) {
      }

      swap(&v[i], &v[k]);
      reverse(j, last, v);
      return true;
    }

    if (i == first) {
      reverse(first, last, v);
      return false;
    }
  }
}

