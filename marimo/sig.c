#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <sys/ipc.h>
#include <sys/shm.h>
#include <math.h>
#include <x86intrin.h> // SIMD命令を使用するためのヘッダファイル


#define M 128
#define K 64
#define N 128
#define E 7

// #include "aaa.c"

#define MATRIX_SIZE 128*64
#define SHM_KEY 1234

short t1[3][6]={{1,0,0,11,6,7},{0,1,0,15,8,14},{0,0,1,2,14,8}};
short t2[6][3]={{11,6,7},{15,8,14},{2,14,8},{1,0,0},{0,1,0},{0,0,1}};

/*
 * S-box transformation table
 */
static const unsigned char s_box[256] = {
    // 0     1     2     3     4     5     6     7     8     9     a     b     c     d     e     f
    0x63, 0x7c, 0x77, 0x7b, 0xf2, 0x6b, 0x6f, 0xc5, 0x30, 0x01, 0x67, 0x2b, 0xfe, 0xd7, 0xab, 0x76,  // 0
    0xca, 0x82, 0xc9, 0x7d, 0xfa, 0x59, 0x47, 0xf0, 0xad, 0xd4, 0xa2, 0xaf, 0x9c, 0xa4, 0x72, 0xc0,  // 1
    0xb7, 0xfd, 0x93, 0x26, 0x36, 0x3f, 0xf7, 0xcc, 0x34, 0xa5, 0xe5, 0xf1, 0x71, 0xd8, 0x31, 0x15,  // 2
    0x04, 0xc7, 0x23, 0xc3, 0x18, 0x96, 0x05, 0x9a, 0x07, 0x12, 0x80, 0xe2, 0xeb, 0x27, 0xb2, 0x75,  // 3
    0x09, 0x83, 0x2c, 0x1a, 0x1b, 0x6e, 0x5a, 0xa0, 0x52, 0x3b, 0xd6, 0xb3, 0x29, 0xe3, 0x2f, 0x84,  // 4
    0x53, 0xd1, 0x00, 0xed, 0x20, 0xfc, 0xb1, 0x5b, 0x6a, 0xcb, 0xbe, 0x39, 0x4a, 0x4c, 0x58, 0xcf,  // 5
    0xd0, 0xef, 0xaa, 0xfb, 0x43, 0x4d, 0x33, 0x85, 0x45, 0xf9, 0x02, 0x7f, 0x50, 0x3c, 0x9f, 0xa8,  // 6
    0x51, 0xa3, 0x40, 0x8f, 0x92, 0x9d, 0x38, 0xf5, 0xbc, 0xb6, 0xda, 0x21, 0x10, 0xff, 0xf3, 0xd2,  // 7
    0xcd, 0x0c, 0x13, 0xec, 0x5f, 0x97, 0x44, 0x17, 0xc4, 0xa7, 0x7e, 0x3d, 0x64, 0x5d, 0x19, 0x73,  // 8
    0x60, 0x81, 0x4f, 0xdc, 0x22, 0x2a, 0x90, 0x88, 0x46, 0xee, 0xb8, 0x14, 0xde, 0x5e, 0x0b, 0xdb,  // 9
    0xe0, 0x32, 0x3a, 0x0a, 0x49, 0x06, 0x24, 0x5c, 0xc2, 0xd3, 0xac, 0x62, 0x91, 0x95, 0xe4, 0x79,  // a
    0xe7, 0xc8, 0x37, 0x6d, 0x8d, 0xd5, 0x4e, 0xa9, 0x6c, 0x56, 0xf4, 0xea, 0x65, 0x7a, 0xae, 0x08,  // b
    0xba, 0x78, 0x25, 0x2e, 0x1c, 0xa6, 0xb4, 0xc6, 0xe8, 0xdd, 0x74, 0x1f, 0x4b, 0xbd, 0x8b, 0x8a,  // c
    0x70, 0x3e, 0xb5, 0x66, 0x48, 0x03, 0xf6, 0x0e, 0x61, 0x35, 0x57, 0xb9, 0x86, 0xc1, 0x1d, 0x9e,  // d
    0xe1, 0xf8, 0x98, 0x11, 0x69, 0xd9, 0x8e, 0x94, 0x9b, 0x1e, 0x87, 0xe9, 0xce, 0x55, 0x28, 0xdf,  // e
    0x8c, 0xa1, 0x89, 0x0d, 0xbf, 0xe6, 0x42, 0x68, 0x41, 0x99, 0x2d, 0x0f, 0xb0, 0x54, 0xbb, 0x16}; // f

/*
 * Inverse S-box transformation table
 */
static const unsigned char inv_s_box[256] = {
    // 0     1     2     3     4     5     6     7     8     9     a     b     c     d     e     f
    0x52, 0x09, 0x6a, 0xd5, 0x30, 0x36, 0xa5, 0x38, 0xbf, 0x40, 0xa3, 0x9e, 0x81, 0xf3, 0xd7, 0xfb,  // 0
    0x7c, 0xe3, 0x39, 0x82, 0x9b, 0x2f, 0xff, 0x87, 0x34, 0x8e, 0x43, 0x44, 0xc4, 0xde, 0xe9, 0xcb,  // 1
    0x54, 0x7b, 0x94, 0x32, 0xa6, 0xc2, 0x23, 0x3d, 0xee, 0x4c, 0x95, 0x0b, 0x42, 0xfa, 0xc3, 0x4e,  // 2
    0x08, 0x2e, 0xa1, 0x66, 0x28, 0xd9, 0x24, 0xb2, 0x76, 0x5b, 0xa2, 0x49, 0x6d, 0x8b, 0xd1, 0x25,  // 3
    0x72, 0xf8, 0xf6, 0x64, 0x86, 0x68, 0x98, 0x16, 0xd4, 0xa4, 0x5c, 0xcc, 0x5d, 0x65, 0xb6, 0x92,  // 4
    0x6c, 0x70, 0x48, 0x50, 0xfd, 0xed, 0xb9, 0xda, 0x5e, 0x15, 0x46, 0x57, 0xa7, 0x8d, 0x9d, 0x84,  // 5
    0x90, 0xd8, 0xab, 0x00, 0x8c, 0xbc, 0xd3, 0x0a, 0xf7, 0xe4, 0x58, 0x05, 0xb8, 0xb3, 0x45, 0x06,  // 6
    0xd0, 0x2c, 0x1e, 0x8f, 0xca, 0x3f, 0x0f, 0x02, 0xc1, 0xaf, 0xbd, 0x03, 0x01, 0x13, 0x8a, 0x6b,  // 7
    0x3a, 0x91, 0x11, 0x41, 0x4f, 0x67, 0xdc, 0xea, 0x97, 0xf2, 0xcf, 0xce, 0xf0, 0xb4, 0xe6, 0x73,  // 8
    0x96, 0xac, 0x74, 0x22, 0xe7, 0xad, 0x35, 0x85, 0xe2, 0xf9, 0x37, 0xe8, 0x1c, 0x75, 0xdf, 0x6e,  // 9
    0x47, 0xf1, 0x1a, 0x71, 0x1d, 0x29, 0xc5, 0x89, 0x6f, 0xb7, 0x62, 0x0e, 0xaa, 0x18, 0xbe, 0x1b,  // a
    0xfc, 0x56, 0x3e, 0x4b, 0xc6, 0xd2, 0x79, 0x20, 0x9a, 0xdb, 0xc0, 0xfe, 0x78, 0xcd, 0x5a, 0xf4,  // b
    0x1f, 0xdd, 0xa8, 0x33, 0x88, 0x07, 0xc7, 0x31, 0xb1, 0x12, 0x10, 0x59, 0x27, 0x80, 0xec, 0x5f,  // c
    0x60, 0x51, 0x7f, 0xa9, 0x19, 0xb5, 0x4a, 0x0d, 0x2d, 0xe5, 0x7a, 0x9f, 0x93, 0xc9, 0x9c, 0xef,  // d
    0xa0, 0xe0, 0x3b, 0x4d, 0xae, 0x2a, 0xf5, 0xb0, 0xc8, 0xeb, 0xbb, 0x3c, 0x83, 0x53, 0x99, 0x61,  // e
    0x17, 0x2b, 0x04, 0x7e, 0xba, 0x77, 0xd6, 0x26, 0xe1, 0x69, 0x14, 0x63, 0x55, 0x21, 0x0c, 0x7d}; // f

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

int reg(short s[N][N])
{

    int count = 0;
    for (int i = 0; i < N; i++)
    {
        if (s[i][i] == 1)
            count++;
    }
    if (count == N)
        return 1;

    return 0;
}
/*
int oinv2(unsigned short b)
{
    int i;

    if (b == 0)
        return 0;

    return (N - fg[b]) % (N - 1) + 1;
}
*/
void T2(short a[N][N],short b[N][N]){
  int i,j,k;

      for(i=0;i<K;i++){
      for(j=0;j<N;j++)
	printf("t%d ",a[i][j]);
      printf("\n");
    }
    //exit(1);
    for(i=K;i<N;i++)
    b[i][i]=1;
    for(i=0;i<K;i++){
      for(j=K;j<N;j++){
	b[j-K][i]=a[i][j];
      }
    }
      for(i=0;i<N;i++){
      for(j=0;j<N;j++)
	printf("%2d ",b[j][i]);
      printf("\n");
    }
    printf("\n");
    exit(1);
}

void make_mat(short a[N][N]){
  int i,j,k;

  for(i=0;i<K;i++){
    a[i][i]=1;
    a[i+K][i]=1;
    for(j=K;j<N;j++){
      a[i][j]=rand()%N;
      a[i+K][j]=rand()%N;
    }
    }
  for(i=0;i<K;i++){
    for(j=0;j<N;j++)
      printf("%d,",a[i][j]);
    printf("\n");
  }
  //exit(1);
}

// #include "gaus.c"

// 白魔法陣 H
short aa2[K][N] = {
    {0, 11, 11, 15, 5, 2, 15, 1, 9, 1, 6, 4, 7, 14, 4, 4},
    {5, 11, 12, 5, 12, 5, 11, 6, 14, 2, 14, 7, 3, 13, 12, 6},
    {1, 15, 12, 6, 13, 10, 10, 15, 8, 11, 0, 13, 1, 5, 15, 12},
    {9, 5, 4, 3, 6, 11, 11, 3, 7, 6, 10, 13, 7, 13, 5, 6},
    {8, 0, 3, 8, 12, 4, 0, 0, 8, 4, 4, 5, 12, 4, 5, 3},
    {0, 15, 11, 1, 15, 14, 12, 13, 4, 3, 10, 10, 3, 12, 2, 1},
    {0, 15, 12, 14, 7, 8, 2, 8, 9, 6, 1, 4, 13, 5, 6, 8},
    {8, 15, 6, 11, 5, 7, 14, 5, 14, 12, 12, 7, 11, 7, 10, 12}};
    /*
    {0, 11, 11, 15, 5, 2, 15, 1, 9, 1, 6, 4, 7, 14, 4, 4},
    {5, 11, 12, 5, 12, 5, 11, 6, 14, 2, 14, 7, 3, 13, 12, 6},
    {1, 15, 12, 6, 13, 10, 10, 15, 8, 11, 0, 13, 1, 5, 15, 12},
    {9, 5, 4, 3, 6, 11, 11, 3, 7, 6, 10, 13, 7, 13, 5, 6},
    {8, 0, 3, 8, 12, 4, 0, 0, 8, 4, 4, 5, 12, 4, 5, 3},
    {0, 15, 11, 1, 15, 14, 12, 13, 4, 3, 10, 10, 3, 12, 2, 1},
    {0, 15, 12, 14, 7, 8, 2, 8, 9, 6, 1, 4, 13, 5, 6, 8},
    {8, 15, 6, 11, 5, 7, 14, 5, 14, 12, 12, 7, 11, 7, 10, 12}};
    */
short a[N][N] = {
    {2, 2, 5, 3, 1, 3, 3, 11, 15, 8, 5, 10, 6, 4, 9, 2},
    {2, 0, 15, 6, 5, 12, 9, 15, 14, 12, 1, 8, 8, 5, 14, 15},
    {0, 0, 13, 10, 7, 9, 3, 11, 9, 6, 10, 7, 10, 12, 11, 10},
    {9, 9, 5, 9, 2, 15, 14, 2, 3, 1, 4, 15, 9, 2, 7, 5},
    {5, 12, 11, 6, 3, 8, 12, 7, 8, 14, 15, 1, 9, 3, 0, 3},
    {5, 9, 4, 4, 7, 0, 0, 15, 15, 0, 13, 9, 9, 14, 12, 7},
    {13, 4, 14, 14, 0, 2, 2, 8, 11, 1, 13, 0, 2, 3, 10, 6},
    {4, 0, 5, 7, 2, 12, 10, 8, 10, 10, 1, 7, 7, 12, 7, 0},
    {6, 13, 12, 15, 6, 3, 9, 11, 5, 4, 5, 12, 4, 10, 7, 8},
    {4, 10, 3, 8, 15, 13, 13, 7, 8, 11, 2, 5, 12, 7, 7, 0},
    {0, 10, 6, 11, 9, 12, 8, 6, 12, 12, 7, 1, 15, 5, 12, 0},
    {6, 7, 0, 1, 4, 10, 12, 10, 15, 2, 6, 7, 12, 6, 3, 8},
    {5, 5, 10, 2, 6, 15, 10, 4, 15, 7, 7, 14, 9, 7, 12, 13},
    {11, 4, 2, 15, 0, 10, 11, 5, 6, 9, 4, 13, 5, 11, 11, 1},
    {0, 4, 4, 2, 0, 4, 15, 8, 5, 13, 14, 6, 9, 6, 8, 15},
    {0, 4, 8, 6, 0, 7, 4, 13, 14, 15, 6, 15, 6, 8, 9, 10}

};

short r2[N][N] = {
    {1, 4, 11, 11, 10, 1, 4, 10, 6, 11, 10, 9, 5, 15, 11, 12},
    {0, 4, 5, 14, 14, 5, 11, 3, 5, 12, 8, 12, 9, 7, 8, 8},
    {2, 12, 15, 4, 10, 0, 7, 14, 2, 3, 8, 12, 12, 8, 12, 10},
    {6, 7, 2, 3, 7, 6, 10, 5, 4, 7, 6, 14, 2, 6, 11, 2},
    {0, 7, 4, 5, 15, 13, 9, 8, 6, 10, 9, 8, 11, 8, 8, 13},
    {15, 14, 11, 12, 5, 3, 10, 1, 1, 6, 9, 9, 11, 8, 10, 15},
    {0, 14, 5, 7, 7, 15, 9, 7, 8, 3, 5, 12, 13, 2, 6, 10},
    {6, 5, 5, 6, 6, 0, 8, 15, 11, 7, 13, 14, 14, 13, 13, 2},
    {1, 4, 11, 11, 10, 1, 4, 10, 6, 11, 10, 9, 5, 15, 11, 12},
    {0, 4, 5, 14, 14, 5, 11, 3, 5, 12, 8, 12, 9, 7, 8, 8},
    {2, 12, 15, 4, 10, 0, 7, 14, 2, 3, 8, 12, 12, 8, 12, 10},
    {6, 7, 2, 3, 7, 6, 10, 5, 4, 7, 6, 14, 2, 6, 11, 2},
    {0, 7, 4, 5, 15, 13, 9, 8, 6, 10, 9, 8, 11, 8, 8, 13},
    {15, 14, 11, 12, 5, 3, 10, 1, 1, 6, 9, 9, 11, 8, 10, 15},
    {0, 14, 5, 7, 7, 15, 9, 7, 8, 3, 5, 12, 13, 2, 6, 10},
    {6, 5, 5, 6, 6, 0, 8, 15, 11, 7, 13, 14, 14, 13, 13, 2}};

// 黒魔法陣 G
short g2[8][16] = {
    {1, 0, 0, 0, 0, 0, 0, 0, 2, 9, 9, 1, 6, 2, 12, 8},
    {0, 1, 0, 0, 0, 0, 0, 0, 6, 7, 2, 13, 13, 5, 6, 3},
    {0, 0, 1, 0, 0, 0, 0, 0, 4, 13, 15, 1, 11, 10, 8, 10},
    {0, 0, 0, 1, 0, 0, 0, 0, 1, 11, 15, 11, 4, 12, 11, 15},
    {0, 0, 0, 0, 1, 0, 0, 0, 15, 9, 14, 10, 11, 9, 11, 11},
    {0, 0, 0, 0, 0, 1, 0, 0, 10, 7, 9, 9, 1, 9, 13, 10},
    {0, 0, 0, 0, 0, 0, 1, 0, 14, 3, 6, 15, 6, 7, 4, 3},
    {0, 0, 0, 0, 0, 0, 0, 1, 13, 5, 2, 11, 7, 14, 10, 8}};
    /*
    {1, 0, 0, 0, 0, 0, 0, 0, 10, 3, 8, 11, 15, 13, 9, 6},
    {0, 1, 0, 0, 0, 0, 0, 0, 12, 5, 13, 13, 13, 13, 12, 9},
    {0, 0, 1, 0, 0, 0, 0, 0, 3, 15, 7, 2, 5, 1, 6, 6},
    {0, 0, 0, 1, 0, 0, 0, 0, 15, 2, 6, 8, 13, 11, 3, 8},
    {0, 0, 0, 0, 1, 0, 0, 0, 12, 6, 2, 4, 1, 15, 2, 7},
    {0, 0, 0, 0, 0, 1, 0, 0, 13, 14, 9, 4, 8, 11, 4, 11},
    {0, 0, 0, 0, 0, 0, 1, 0, 7, 1, 1, 7, 8, 13, 2, 15},
    {0, 0, 0, 0, 0, 0, 0, 1, 4, 10, 2, 7, 13, 3, 13, 2}};
    */
short ggg[16][16] = {
    {1, 0, 0, 0, 0, 0, 0, 0, 3, 13, 8, 2, 6, 7, 13, 1},
    {0, 1, 0, 0, 0, 0, 0, 0, 4, 8, 9, 5, 11, 14, 12, 7},
    {0, 0, 1, 0, 0, 0, 0, 0, 1, 11, 14, 14, 4, 13, 6, 2},
    {0, 0, 0, 1, 0, 0, 0, 0, 3, 11, 4, 10, 8, 2, 1, 7},
    {0, 0, 0, 0, 1, 0, 0, 0, 12, 3, 12, 2, 13, 9, 6, 7},
    {0, 0, 0, 0, 0, 1, 0, 0, 2, 13, 11, 3, 3, 15, 10, 7},
    {0, 0, 0, 0, 0, 0, 1, 0, 3, 5, 9, 8, 7, 10, 9, 12},
    {0, 0, 0, 0, 0, 0, 0, 1, 1, 12, 3, 11, 10, 3, 1, 9},
    {1, 0, 0, 0, 0, 0, 0, 0, 3, 13, 8, 2, 6, 7, 13, 1},
    {0, 1, 0, 0, 0, 0, 0, 0, 4, 8, 9, 5, 11, 14, 12, 7},
    {0, 0, 1, 0, 0, 0, 0, 0, 1, 11, 14, 14, 4, 13, 6, 2},
    {0, 0, 0, 1, 0, 0, 0, 0, 3, 11, 4, 10, 8, 2, 1, 7},
    {0, 0, 0, 0, 1, 0, 0, 0, 12, 3, 12, 2, 13, 9, 6, 7},
    {0, 0, 0, 0, 0, 1, 0, 0, 2, 13, 11, 3, 3, 15, 10, 7},
    {0, 0, 0, 0, 0, 0, 1, 0, 3, 5, 9, 8, 7, 10, 9, 12},
    {0, 0, 0, 0, 0, 0, 0, 1, 1, 12, 3, 11, 10, 3, 1, 9}};

int mlt2(int x, int y)
{

    if (x == 0 || y == 0)
        return 0;

    return ((x + y - 2) % (M - 1)) + 1;
}


void mat_print(short A[N][N])
{
    int i, j;
    for (i = 0; i < N; i++)
    {
        for (j = 0; j < N; j++)
            printf("%d,", A[i][j]);
        printf("\n");
    }
    printf("\n");
}



// 行列を表示する関数
void print_matrix(short A[MATRIX_SIZE][MATRIX_SIZE])
{
    for (int i = 0; i < MATRIX_SIZE; i++)
    {
        for (int j = 0; j < MATRIX_SIZE; j++)
        {
            printf("%d ", A[i][j]);
        }
        printf("\n");
    }
}



// 行列の逆行列を計算する関数
void inverseMatrix(short A[MATRIX_SIZE][MATRIX_SIZE], short A_inv[MATRIX_SIZE][MATRIX_SIZE], int start_row, int end_row)
{
    int i, j, k;
    short temp;

    // 単位行列を初期化
    for (i = start_row; i < end_row; i++)
    {
        for (j = 0; j < MATRIX_SIZE; j++)
        {
            A_inv[i][j] = (i == j) ? 1.0 : 0.0;
        }
    }

    // ガウス・ジョルダン法による逆行列の計算
    for (k = start_row; k < end_row; k++)
    {
        temp = A[k][k];
        for (j = 0; j < MATRIX_SIZE; j++)
        {
            A[k][j] /= temp;
            A_inv[k][j] /= temp;
        }
        for (i = start_row; i < end_row; i++)
        {
            if (i != k)
            {
                temp = A[i][k];
                for (j = 0; j < MATRIX_SIZE; j++)
                {
                    A[i][j] -= A[k][j] * temp;
                    A_inv[i][j] -= A_inv[k][j] * temp;
                }
            }
        }
    }
}


void Ta(short H[K][K], short G[K][N])
{
    int i, j;

    for (i = 0; i < K; i++)
        G[i][i] = 1;
    for (i = 0; i < K; i++)
    {
        for (j = K; j < N; j++)
            G[j - K][i + K] = H[i][j - K];
    }
}




long long sum(short a[N]){
    long long u=0;
    for(int i=0;i<N;i++)
    if(a[i]+a[(i+1)%N]!=0)
        u+=(a[i]^a[(i+1)%N]);

    return u;
}

int main()
{
    short A[MATRIX_SIZE][MATRIX_SIZE];
    short A_inv[MATRIX_SIZE][MATRIX_SIZE];
    short C[MATRIX_SIZE][MATRIX_SIZE];
    short AA[MATRIX_SIZE][MATRIX_SIZE];

    short rr[N][N] = {0};
    short ss[N][N] = {0};
    short tt[N][N] = {0};
    short ga[N][N] = {0};
    short z[N]={0};
    short v[N]={0};
    short ra[N]={0};
    short w[N]={0};

    int i, j, k;
    short c[N][N] = {0};
    short d[K][K] = {0};
    short P[N] = {0},Pa[N]={0};
    short bb[N][N] = {0};
    short cc[K][N] = {0};
    short cd[N][N] = {0};
    short b[N][N] = {0};
    short a2[N][N] = {0};
    short a3[N][N] = {0};

    memcpy(a2, a, sizeof(a2));
    memcpy(a3, a, sizeof(a2));

    // vis();
    // exit(1);
    srand(clock());

    for (i = 0; i < N; i++){
        P[i] = i;
        Pa[i]=i;
    }
    random_shuffle(P, N);
    random_shuffle(Pa,N);
    short inv_P[N];
    for (i = 0; i < N; i++)
    {
        printf("%d,", P[i]);
        inv_P[P[i]] = i;
    }
    short inv_Pa[N]={0};
    for(i=0;i<N;i++)
    inv_Pa[Pa[i]]=i;

    printf("\n");
    for(i=0;i<N;i++){
    ra[i]=rand()%4096;
    w[i]=rand()%4096;
    if(rand()%2==1)
    ra[i]= -1*ra[i];
    if(rand()%2==1)
    v[i]= -1*w[inv_Pa[i]];
    }
    
    short p1[N]={0},p2[N]={0},p3[N]={0};
    int ca=517;
    for(i=0;i<N;i++){
    p1[i]=ra[i]+v[i];

    p2[i]=ra[P[i]];
    p3[i]=ra[P[i]]+v[P[i]];
    //z[i]=ca*v[Pa[i]];
    }
    
    for(i=0;i<N;i++)
    printf("%d %d\n",p1[i],p2[i]);
    long long c0=sum(p1);
    long long c1=sum(p2);
    printf("%lld\n",c0);
    printf("%lld\n",c1);
    
    for(i=0;i<N;i++){
    z[i]=ra[P[i]]+ca*v[Pa[P[i]]];
    //printf("aa%d %d\n",z[i],p3[i]);
    }
    short va[N]={0};
    for(i=0;i<N;i++)
    va[i]=p3[inv_P[i]];
    printf("%lld\n",sum(va));
    short vb[N]={0};
    for(i=0;i<N;i++)
    vb[i]=z[i]-ca*v[Pa[P[i]]];
    printf("%lld\n",sum(vb));


    return 0;
}
