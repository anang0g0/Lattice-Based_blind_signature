#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <inttypes.h>
#include <immintrin.h>

unsigned char p[32], q[32],r[32],inv_r[32],inv_p[32],gf[256],fg[256],kkk[32]={0};

static const uint8_t s_box[256]={198,199,70,57,92,213,163,101,243,107,133,248,150,155,108,181,3,216,45,145,16,53,244,135,83,192,137,230,128,240,99,222,167,22,95,59,77,30,254,73,228,161,249,242,174,20,11,206,50,55,235,34,74,219,109,1,166,134,152,239,65,232,114,189,160,97,201,149,104,229,184,153,171,113,165,121,217,82,157,33,118,39,141,88,116,31,131,93,76,215,210,86,203,119,170,56,84,245,226,117,183,48,140,12,6,209,196,190,0,236,188,204,32,127,139,251,18,29,129,144,241,9,44,96,25,172,15,63,13,225,90,194,7,125,200,247,182,193,246,110,185,69,146,176,250,75,130,94,187,234,238,227,223,36,178,191,164,162,23,214,47,173,58,103,124,2,197,85,52,64,37,21,61,168,115,147,43,154,158,233,40,19,132,72,28,138,175,100,122,123,35,78,159,156,46,186,91,10,180,67,120,98,79,89,252,169,102,14,17,5,179,41,221,237,148,27,60,224,26,211,143,106,177,112,151,42,195,66,81,212,111,231,255,54,62,80,38,202,126,4,24,220,208,49,205,71,218,68,8,105,87,136,253,207,142,51};
static const uint8_t inv_s_box[256]={108,55,165,16,239,209,104,132,248,121,197,46,103,128,207,126,20,208,116,181,45,171,33,158,240,124,218,215,184,117,37,85,112,79,51,190,153,170,236,81,180,211,225,176,122,18,194,160,101,243,48,255,168,21,233,49,95,3,162,35,216,172,234,127,169,60,227,199,247,141,2,245,183,39,52,145,88,36,191,202,235,228,77,24,96,167,91,250,83,203,130,196,4,87,147,34,123,65,201,30,187,7,206,163,68,249,221,9,14,54,139,230,223,73,62,174,84,99,80,93,200,75,188,189,164,133,238,113,28,118,146,86,182,10,57,23,251,26,185,114,102,82,254,220,119,19,142,175,214,67,12,224,58,71,177,13,193,78,178,192,64,41,157,6,156,74,56,32,173,205,94,72,125,161,44,186,143,222,154,210,198,15,136,100,70,140,195,148,110,63,107,155,25,137,131,226,106,166,0,1,134,66,237,92,111,244,47,253,242,105,90,219,229,5,159,89,17,76,246,53,241,212,31,152,217,129,98,151,40,69,27,231,61,179,149,50,109,213,150,59,29,120,43,8,22,97,138,135,11,42,144,115,204,252,38,232};

int Nr;
int Nk;
int Round;

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


uint16_t inv(uint16_t a,uint16_t n)
{
  uint16_t d = n;
  uint16_t x = 0;
  uint16_t s = 1,q,r,t;
  while (a != 0){
    q = d / a;
    r = d % a;
    d = a;
    a = r;
    t = x - q * s;
    x = s;
    s = t;
  }
  //gcd = d  # $\gcd(a, n)$ 

  return ((x + n) % (n / d));
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
uint8_t gmult(uint8_t a, uint8_t b)
{

    uint8_t p = 0, i = 0, hbs = 0;

    for (i = 0; i < 8; i++)
    {
        if (b & 1)
        {
            p ^= a;
        }

        hbs = a & 0x80;
        a <<= 1;
        if (hbs)
            a ^= 0x1b; // 0000 0001 0001 1011
        b >>= 1;
    }

    return (uint8_t)p;
}

void ha(__uint128_t a,__uint128_t b){
    __uint128_t p = {0}, hbs = {0},out=1;

    for (int i = 0; i < 8; i++)
    {
        if (b & out)
        {
            p ^= a;
        }

        hbs = a & 0x80;
        a <<= 1;
        if (hbs)
            a ^= 0xf3; // 0000 0000 1111 0011
        b >>= 1;
    }

}


/*
 * Generates the round constant Rcon[i]
 */
uint8_t R[] = {0x02, 0x00, 0x00, 0x00};
 
uint8_t * Rcon(uint8_t i) {
	
	if (i == 1) {
		R[0] = 0x01; // x^(1-1) = x^0 = 1
	} else if (i > 1) {
		R[0] = 0x02;
		i--;
		while (i > 1) {
			R[0] = gmult(R[0], 0x02);
			i--;
		}
	}
	
	return R;
}

#define Nb 8
/*
 * Transformation in the Cipher and Inverse Cipher in which a Round 
 * Key is added to the State using an XOR operation. The length of a 
 * Round Key equals the size of the State (i.e., for Nb = 4, the Round 
 * Key length equals 128 bits/16 bytes).
 */
void add_round_key(uint8_t *state, uint8_t *w, uint8_t r) {
	
	uint8_t c;
	
	for (c = 0; c < Nb; c++) {
		state[Nb*0+c] = state[Nb*0+c]^w[4*Nb*r+4*c+0];   //debug, so it works for Nb !=4 
		state[Nb*1+c] = state[Nb*1+c]^w[4*Nb*r+4*c+1];
		state[Nb*2+c] = state[Nb*2+c]^w[4*Nb*r+4*c+2];
		state[Nb*3+c] = state[Nb*3+c]^w[4*Nb*r+4*c+3];	
	}
}

/*
 * Multiplication of 4 byte words
 * m(x) = x4+1
 */
void coef_mult(uint8_t *a, uint8_t *b, uint8_t *d) {

	d[0] = gmult(a[0],b[0])^gmult(a[3],b[1])^gmult(a[2],b[2])^gmult(a[1],b[3]);
	d[1] = gmult(a[1],b[0])^gmult(a[0],b[1])^gmult(a[3],b[2])^gmult(a[2],b[3]);
	d[2] = gmult(a[2],b[0])^gmult(a[1],b[1])^gmult(a[0],b[2])^gmult(a[3],b[3]);
	d[3] = gmult(a[3],b[0])^gmult(a[2],b[1])^gmult(a[1],b[2])^gmult(a[0],b[3]);
}


/*
 * Function used in the Key Expansion routine that takes a four-byte 
 * word and performs a cyclic permutation. 
 */
void rot_word(uint8_t *w) {

	uint8_t tmp;
	uint8_t i;

	tmp = w[0];

	for (i = 0; i < 3; i++) {
		w[i] = w[i+1];
	}

	w[3] = tmp;
}

/*
 * Addition of 4 byte words
 * m(x) = x4+1
 */
void coef_add(uint8_t a[], uint8_t b[], uint8_t d[]) {

	d[0] = a[0]^b[0];
	d[1] = a[1]^b[1];
	d[2] = a[2]^b[2];
	d[3] = a[3]^b[3];
}

/*
 * Function used in the Key Expansion routine that takes a four-byte 
 * input word and applies an S-box to each of the four bytes to 
 * produce an output word.
 */
void sub_word(uint8_t *w) {

	uint8_t i;

	for (i = 0; i < 4; i++) {
		w[i] = s_box[w[i]];
	}
}

/*
 * Key Expansion
 */
void aes_key_expansion(uint8_t *key, uint8_t *w) {

	uint8_t tmp[4];
	uint8_t i;
	uint8_t len = Nb*(Nr+1);

	for (i = 0; i < Nk; i++) {
		w[4*i+0] = key[4*i+0];
		w[4*i+1] = key[4*i+1];
		w[4*i+2] = key[4*i+2];
		w[4*i+3] = key[4*i+3];
	}

	for (i = Nk; i < len; i++) {
		tmp[0] = w[4*(i-1)+0];
		tmp[1] = w[4*(i-1)+1];
		tmp[2] = w[4*(i-1)+2];
		tmp[3] = w[4*(i-1)+3];

		if (i%Nk == 0) {

			rot_word(tmp);
			sub_word(tmp);
			coef_add(tmp, Rcon(i/Nk), tmp);

		} else if (Nk > 6 && i%Nk == 4) {

			sub_word(tmp);

		}

		w[4*i+0] = w[4*(i-Nk)+0]^tmp[0];
		w[4*i+1] = w[4*(i-Nk)+1]^tmp[1];
		w[4*i+2] = w[4*(i-Nk)+2]^tmp[2];
		w[4*i+3] = w[4*(i-Nk)+3]^tmp[3];
	}
}



void rounder(){
    for(int i=0;i<32;i++)
    q[i]=p[r[inv_p[i]]];
	memcpy(r,q,32);
	for(int i=0;i<32;i++)
	inv_r[r[i]]=i;
}

void reverse(){
    for(int i=0;i<32;i++)
    q[i]=inv_p[r[p[i]]];
	memcpy(r,q,32);
	for(int i=0;i<32;i++)
	inv_r[r[i]]=i;
}

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


int mlt(int x, int y)
{

  if (x == 0 || y == 0)
    return 0;

  return ((x + y - 2) % (256 - 1)) + 1;
}


void add(uint8_t *m,uint8_t *k){
    int i,n;
    for (i = 0; i < 32; i++){
        n = (m[i] + k[r[i]]) % 256;
		//m[i]=s_box[((n % 16) + (n >> 4) * 16)];
        m[i]=n;
    }
}

void sub(uint8_t *c,uint8_t *k){
    int i,n;
    for(i=0;i<32;i++){
        c[i]=c[i];
        //c[i]=inv_s_box[((n % 16) + (n >> 4) * 16)];
        c[i] = (256+ c[i] - k[r[i]]) % 256;
    }

}


void perm(uint8_t *m,uint8_t *r){
	uint8_t u[32];
	for(int i=0;i<32;i++)
	u[i]=m[r[i]];
	memcpy(m,u,32);
}

void sche(uint8_t *o){
	for(int i=0;i<32;i++)
	o[i]^=o[r[i]];
}

void euc(uint8_t *m,uint8_t *k){
int i,j;
uint8_t s[16],tmp[16],u[32];

for(i=0;i<16;i++){
tmp[i]=m[i];
s[i]=s_box[m[i]^k[i]]^m[i+16];
m[i+16]=tmp[i];
m[i]=s[i];
}
}

void uec(uint8_t *m,uint8_t *k){
int i,j;
uint8_t s[16],tmp[16],u[32]={0};

for(i=0;i<16;i++){
tmp[i]=m[i];
m[i]^=s_box[m[i+16]];
u[i+16]=m[i];
m[i]=inv_s_box[tmp[i]^u[i+16]];
m[i+16]=u[i+16];
}
}


void enc(uint8_t *m, __uint8_t *k)
{
    int i;
    unsigned char n=0,ff[32]={0};
	

    //rounder();
	//add_round_key(m,k,Round);
	//sche(k);
		//memcpy(k,ff,32);
    for (i = 0; i < 4; i++)
	{
		unsigned char tmp[32]={0};

		for(int j=0;j<8;j++){
		//k[i]^=k[inv_r[i*8+j]];
		m[i*8+j]^=k[inv_r[(i*8+j)]];
		n= (m[(i*8+j)%32] + (k[r[(i*8+j)%32]])) % 256;
		m[(i*8+j)]=s_box[((n % 16) + (n >> 4) * 16)];
		}

    }
	uint8_t tmp[32];
	for(int i=0;i<32;i++)
	tmp[i]=m[r[i]];
	//memcpy(m,tmp,32);
	//for(int l=0;l<32;l++)
	//k[l]^=k[r[l]];

    //shift_rows(m);
    //mix_columns(m);
}

void dec(uint8_t *c, uint8_t *k)
{
    int i;
    unsigned char n=0,ff[32]={0};

    //inv_mix_columns(c);
    //inv_shift_rows(c);

	uint8_t tmp[32];
	for(i=0;i<32;i++)
	tmp[i]=c[inv_r[i]];
	//memcpy(c,tmp,32);
	
	//add_round_key(c,k,Round);
    for (i = 0; i < 4; i++){
		unsigned char tmp[32]={0};
		
		for(int j=0;j<8;j++){
		n=c[i*8+j];
	    c[i*8+j]=inv_s_box[((n % 16) + (n >> 4) * 16)];
		c[(i*8+j)%32] = (256 + c[(i*8+j)%32] - (k[r[(i*8+j)%32]])) % 256;
		//k[i]^=k[inv_r[i*8+j]];
		c[i*8+j]^=k[inv_r[i*8+j]];
		}
    }

	//for(int l=0;l<32;l++)
	//k[l]^=k[inv_r[l]];
		//
		//memcpy(k,ff,32);
	//sche(k);
	//reverse();
}


/*
 * Performs the AES cipher operation
 
void aes_cipher(uint8_t *in, uint8_t *out, uint8_t *w) {

	uint8_t state[4*Nb];
	uint8_t r, i, j;

	for (i = 0; i < 4; i++) {
		for (j = 0; j < Nb; j++) {
			state[Nb*i+j] = in[i+4*j];
		}
	}

	add_round_key(state, w, 0);

	for (r = 1; r < Nr; r++) {
		sub_bytes(state);
		shift_rows(state);
		mix_columns(state);
		add_round_key(state, w, r);
	}

	sub_bytes(state);
	shift_rows(state);
	add_round_key(state, w, Nr);

	for (i = 0; i < 4; i++) {
		for (j = 0; j < Nb; j++) {
			out[i+4*j] = state[Nb*i+j];
		}
	}
}
*/
/*
 * Performs the AES inverse cipher operation
 
void aes_inv_cipher(uint8_t *in, uint8_t *out, uint8_t *w) {

	uint8_t state[4*Nb];
	uint8_t r, i, j;

	for (i = 0; i < 4; i++) {
		for (j = 0; j < Nb; j++) {
			state[Nb*i+j] = in[i+4*j];
		}
	}

	add_round_key(state, w, Nr);

	for (r = Nr-1; r >= 1; r--) {
		inv_shift_rows(state);
		inv_sub_bytes(state);
		add_round_key(state, w, r);
		inv_mix_columns(state);
	}

	inv_shift_rows(state);
	inv_sub_bytes(state);
	add_round_key(state, w, 0);

	for (i = 0; i < 4; i++) {
		for (j = 0; j < Nb; j++) {
			out[i+4*j] = state[Nb*i+j];
		}
	}
}
*/

/*
 * Initialize AES variables and allocate memory for expanded key
 */
uint8_t *aes_init(size_t key_size) {

        switch (key_size) {
		default:
		case 16: Nk = 4; Nr = 10; break;
		case 24: Nk = 6; Nr = 12; break;
		case 32: Nk = 8; Nr = 14; break;
	}

	return malloc(Nb*(Nr+1)*4);
}

unsigned short mkbox(unsigned short i){
	return (198^(i*((i*i)%257*(i*i)%257)%257))%257;
}

int invb(int i){
	return inv(i,257);
}

uint64_t pd(uint64_t a, uint64_t b,uint64_t d)
{
  uint64_t c = 0,hbs=0,ll=(1<<(8-1)), l=d^(1<<8);

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

void gen_t_box(unsigned *gf,unsigned *fg)
{
  uint32_t i, j;

gf[0]=0;
  for (i = 0; i < 256; i++)
  {
    gf[i]=pd(seki(seki(seki(i,i),seki(i,i)),seki(i,i)),i,0b100011011)^0b11000110; //pmod(pmod(i,i,normal[E]),i,normal[E]);
  }

  //exit(1);
  printf("static const unsigned short gf[256]={\n");
  for (i = 0; i < 256; i++)
  {
    fg[gf[i]]=i;
    printf("%3d,", gf[i]);
  }
  printf("};\n");
  printf("\n");
  printf("static const unsigned short fg[256]={\n");
  for (i = 0; i < 256; i++)
  {
    printf("%3d,", fg[i]);
  }
  printf("};\n");
  printf("\n");
}


int main()
{
    unsigned short i;
    unsigned char m[128];
    unsigned char k[32],ss[32];
	unsigned tt[256],inv_tt[256];
    //unsigned char p[32],inv_p[32];
	unsigned char s[32],inv_s[32];
    for (i = 0; i < 32; i++)
    {
        m[i] = i+1;
        k[i] = rand() % 256;
        p[i]=i;
        r[i]=i;
		ss[i]=k[i];
		
    }
	uint8_t *w; // expanded key
    uint8_t out[32];

	//exit(1);
	gen_t_box(tt,inv_tt);

	for(i=0;i<256;i++){
	printf("%u,",tt[i]);
	inv_tt[tt[i]]=i;
	}
	printf("\n");
	for(i=0;i<256;i++)
	printf("%u,",inv_tt[i]);
	printf("\n");
	//exit(1);
	
	w = aes_init(32);

	//aes_key_expansion(k, w);
	memcpy(w,k,32);


    random_shuffle(p,32);
    random_shuffle(r,32);
	memcpy(s,r,32);
    for(i=0;i<32;i++){
        inv_p[p[i]]=i;
		inv_r[r[i]]=i;
    }
	memcpy(out,r,32);
	//for(i=0;i<32;i++)
	//printf("%d ",r[i]);
	printf("\n");
	printf("Plaintext message:\n");
	for (i = 0; i < 32; i++) {
		printf("%02x ", m[i]);
	}
	printf("\n");
	uint8_t table[16][32];
	for(i=0;i<16;i++){
	for(int j=0;j<32;j++)
	w[j]^=w[r[j]];
	rounder();
	for(int j=0;j<32;j++)
	table[i][j]=w[j];
	}
	memcpy(r,out,32);
	//aes_cipher(m /* in */, out /* out */, w /* expanded key */);
	for(i=0;i<16;i++){
	for(int l=0;l<32;l++){
	printf("%d ",w[l]);
	//w[l]^=w[r[l]];
	w[l]=table[i][l];
	}
	printf("\n");
	rounder();
	//perm(m,r);
	//enc(m,w);
	add(m,w);
	//euc(m,w);
	for(int l=0;l<32;l++)
	m[l]=s_box[m[l]];
	}

	printf("Ciphered message:\n");
	for (i = 0; i < 32; i++) {
		printf("%02x ", m[i]);
	}
	printf("\n");
	//exit(1);
	
	//memcpy(w,ss,32);
	memcpy(r,out,32);
	for(i=0;i<16;i++)
	rounder();
	for(i=0;i<32;i++)
	printf("%d ",r[i]);
	printf("\n");

	//aes_inv_cipher(out, m, w);
	for(i=0;i<16;i++){
	for(int l=0;l<32;l++){
	printf("%d ",w[l]);
	w[l]=table[15-i][l];
	}
	printf("\n");
	for(int l=0;l<32;l++)
	m[l]=inv_s_box[m[l]];
	//uec(m,w);
	sub(m,w);
	//for(int j=0;j<32;j++)
	//w[j]=table[15-i][j];
	//dec(m,w);
	//perm(m,inv_r);
	//rounder();
	reverse();
	}

	printf("Original message (after inv cipher):\n");
	for (i = 0; i < 32; i++) {
		printf("%02x ", inv_s_box[s_box[m[i]]]);
	}

	printf("\n");

	free(w);

}

/*
	for(int l=0;l<16;l++){
	m[l]=s_box[m[l]];
	//tmp[l]=m[l];
	//m[l]=s_box[m[l]]^m[l+16];
	//m[l+16]=tmp[l];
	}
	}

	printf("Ciphered message:\n");
	for (i = 0; i < 32; i++) {
		printf("%02x ", m[i]);
	}

	printf("\n");
	memcpy(w,ss,32);
	memcpy(r,out,32);
	//for(i=0;i<32;i++)
	//printf("%d ",r[i]);
	//printf("\n");

	//aes_inv_cipher(out, m, w);
	for(i=0;i<16;i++){
	for(int l=0;l<32;l++){
	printf("%d ",w[l]);
	w[l]^=w[r[l]];
	}
	printf(" bbb\n");
	
	for(int l=0;l<16;l++){
	m[l]=inv_s_box[m[l]];
	//tmp[l]=m[l+16];
	//m[l+16]=s_box[m[l+16]]^m[l];
	//m[l]=tmp[l];
	}
*/	
