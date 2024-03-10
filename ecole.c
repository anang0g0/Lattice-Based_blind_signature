/* generate GF(2^n) using irreducible polynomial */
// ゼフ対数表を作るためのプログラム。正規基底を生成します。

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#define O 256
#define E 8

/* generate Galois Field over GF(2^?) */
static const unsigned long long int normal[17] = {
  0b0,
  0b0,
  0b0,
    0b1011,
  //0b11001", /* GF(16) */
    0b10011,
    0b110111,
    0b1100001,
    0b11000001,
 // 0b110101001,
  //  0b100011101, // sage
    0b100011011, // aes
 // 0b100011011",
    0b1100110001,
 // 0b11000010011,
    0b10001101111, // sage1024
    0b110000001101,
    0b1000011101011, // sage 4096
 // 0b1100101000001, /* 4096 */
    0b11011000000001, /* 8192 */
  //  0b10000000011011, /* Classic McEliece */
    0b110000100010001,
    0b1100000000000001,
    0b11010000000010001};

unsigned int gf[O], fg[O];

uint32_t pdiv(uint32_t a, uint32_t b)
{
  uint32_t c = 0,d=a;

  //printf("%b,%b\n", a, b);
  // while(b>0)
  if(a==b){
    printf("eq\n");
    //exit(1);
  return 0;
  }
  c = b;
  while (1)
  {
    int m=a,n=b,count=0,cnt=0;
    while(m>0){
      cnt++;
      m=m>>1;
    }
    while(n>0){
      count++;
      n=n>>1;
    }
    while (cnt > count){
      b = (b << 1);
      count++;
    }
    if(cnt==count && a%2==1){
    //printf("545646 %b %b %b\n",a,b,(b^a));
    a^=(b);
    if(a<c)
    return a&0xff;
    }
    if(a%2==0 && cnt==count)
    a ^= c;
    if (c > a)
      break;

    b = c;
  }

  return a&0xff;
}

uint32_t pd(uint32_t a, uint32_t b,uint32_t d)
{
  uint32_t c = 0,hbs=0,ll=(1<<(E-1)), l=d^(1<<E);

  while (a != 0)
  {
    //printf("b %b %b %b\n",a,b,c);
    if ((a & 1) == 1)
      c ^= b;
    hbs= b&(ll);
    b <<= 1;
    if(hbs)
    b ^=l;
    a >>= 1;
  }

  return c&0xff; //mask
}

uint64_t seki(uint64_t a, uint64_t b)
{
  uint64_t c = 0,hbs=0;

  while (a != 0)
  {
    //printf("b %b %b %b\n",a,b,c);
    if ((a & 1) == 1)
      c ^= b;
    
    b <<= 1;
    
    a >>= 1;
  }

  return c;
}

uint32_t pmod(uint32_t a, uint32_t b, uint32_t c)
{

  //printf("vv %b %b %b\n",seki(a,b),c,pdiv(seki(a,b),c));
  return pdiv(seki(a, b), c);
}

/*
 * Multiplication in GF(2^8)
 * http://en.wikipedia.org/wiki/Finite_field_arithmetic
 * Irreducible polynomial m(x) = x8 + x4 + x3 + x + 1
 *
 * NOTE: This function can be easily replaced with a look up table for a speed
 *       boost, at the expense of an increase in memory size (around 65 KB). See
 *       the aes.h header file to find the macro definition.
 *
 */
uint32_t gmult(uint32_t a, uint32_t b)
{

  uint32_t p = 0, i = 0, hbs = 0;

  for (i = 0; i < 9; i++)
  {
    if (b & 1)
    {
      p ^= a;
    }

    hbs = a & (0x80);
    a <<= 1;
    if (hbs)
      a ^= 0x1b; // 0000 0001 0001 1011
    b >>= 1;
  }

  return (uint8_t)p;
}

void makefg(int n)
{
  int i, j;

  for (i = 0; i < n; i++)
  {
    for (j = 0; j < n; j++)
    {
      if (gf[i] == j)
        fg[j] = i;
    }
  }
  printf("static const unsigned short fg[%d]={", O);
  for (i = 0; i < O; i++)
  {
    if (i < O - 1)
    {
      printf("%d,", fg[i]);
    }
    else
    {
      printf("%d", fg[i]);
    }
  }
  printf("};\n");

  return;
}

void mkgf(int n)
{
  int i, j, bit, count = 0;
  unsigned int pol, N, M, L;

  for (i = 0; i < 13; i++)
    pol = normal[i]; // strtoul(normal[i],(char **)NULL,2);

  /* define pol */
  switch (n)
  {

  case 8:
    pol = normal[3];
    printf("%d\n", n);
    break;

  case 16:
    pol = normal[4];
    printf("%d\n", n);
    break;

  case 32:
    pol = normal[5];
    printf("%d\n", n);
    break;

  case 64:
    pol = normal[6];
    printf("%d\n", n);
    break;

  case 128:
    pol = normal[7];
    printf("%d\n", n);
    break;

  case 256:
    pol = normal[8];
    printf("%d\n", n);
    break;

  case 512:
    pol = normal[9];
    printf("%d\n", n);
    break;

  case 1024:
    pol = normal[10];
    printf("%d\n", n);
    break;

  case 2048:
    pol = normal[11];
    printf("%d\n", n);
    break;

  case 4096:
    pol = normal[12];
    printf("%d\n", n);
    break;

  case 8192:
    pol = normal[13];
    printf("%d\n", n);
    break;

  default: /* 16384 */
    pol = normal[14];
    printf("%d\n", n);
    break;
  }
  L = 1;
  while (pol > L) // 原始多項式の最大次数を計算する。
  {
    L = (L << 1);
    count++;
  }
  L = (L >> 1);
  N = pol ^ L; // 原始多項式の最大次数を消した多項式の残り。

  gf[0] = 0;
  bit = 1;
  for (i = 1; i < L; i++)
  {
    if (bit > L - 1) // もしbitが最大次数に達したら
    {
      bit = bit - L; // bitから最大次数の項 x^n を消す。
      bit = bit ^ N; // 最大次数の項以下の原始多項式を bit に xorする。
    }
    gf[i] = bit;      // 最初は何もしないでx^iを入れていく。
    bit = (bit << 1); // 原始多項式の次数を1上げる。
  }
  printf("static const unsigned short gf[%d]={", O);
  for (i = 0; i < L; i++)
  {
    if (i < O - 1)
    {
      printf("%u,", gf[i]);
    }
    else
    {
      printf("%u", gf[i]);
    }
  }

  printf("};\n");
}

void make()
{
  uint32_t i, j;

  printf("static const unsigned short fg[256]={\n");
  gf[0] = 0;
  printf("  0,");
  for (i = 1; i < O; i++)
  {
    for (j = 0; j < O; j++)
    {
    //
    //
    //if(seki(i,j)==1)
    //if (pd(i, j) == 1)
      //printf("iiiii %b %b %b\n",i,j,pmod(i,j,normal[13]));
      if(pmod(i,j,normal[E])==1)
      {
        gf[i] = j;
        printf("---%3d %3d %b %b %b,\n", i,j,pmod(i,j,normal[E]),pd(i,j,normal[E]),seki(i,j));
      }
    }
  }
  printf("};\n");
  //exit(1);
  for (i = 0; i < O; i++)
  {
    printf("%3d,", gf[i]);
  }
  printf("\n");
}

void t_box()
{
  uint32_t i, j;


  gf[0] = 0;
  for (i = 1; i < O; i++)
  {
    gf[i]=gmult(gmult(i,i),i); //pmod(pmod(i,i,normal[E]),i,normal[E]);
  }

  //exit(1);
  printf("static const unsigned short gf[256]={\n");
  for (i = 0; i < O; i++)
  {
    //fg[gf[i]]=i;
    printf("%3d,", gf[i]);
  }
  printf("};\n");
  printf("\n");
  printf("static const unsigned short fg[256]={\n");
  for (i = 0; i < O; i++)
  {
    //printf("%3d,", fg[i]);
  }
  printf("};\n");
  printf("\n");
}

int main()
{
  int i, j, k;

  printf("%b %b %b %b %b\n",4,5, pmod(4, 0b101,normal[3]), pd(0b100, 0b101,normal[3]),(seki(4,5)^(normal[3]<<1)));
    printf("%b %b %b %b\n",4,7, pmod(4, 0b111,normal[3]), pd(0b100, 0b111,normal[3]));
  //printf("%b %b %b\n", pdiv(0b101111,0b101),pmod(0b101111,0b101,normal[13]),pd(0b101111,2,normal[13]));
  //printf("%b %b\n", pd(6912, 0b10, normal[13])); //pmod(0b100000000,0b10,0b100011011));
  //exit(1);

  make();
  //t_box();
  printf("%d %d\n",gmult(seki(9,9),9),gmult(seki(34,34),34));
  exit(1);

  mkgf(O);
  makefg(O);

  return 0;
}
