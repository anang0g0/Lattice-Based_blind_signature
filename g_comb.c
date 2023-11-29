/*
 * Author: Hiroyuki Chishiro
 * License: 2-Clause BSD
 */
#include <stdio.h>
#include <stdbool.h>

// modified. use ./a.out|grep @
#define N 40

unsigned char xx[256][23] = {0};

void swap(unsigned char *pa, unsigned char *pb)
{
	int tmp;

	tmp = *pa;
	*pa = *pb;
	*pb = tmp;
}

void rotate(size_t first, size_t middle, size_t last, unsigned char v[])
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

bool next_combination(size_t first, size_t last, size_t r, unsigned char v[])
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
	unsigned char v[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193};

	size_t i;
	size_t n = 23; // sizeof(v) / sizeof(v[0]);
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