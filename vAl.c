#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <inttypes.h>
#include <immintrin.h>

unsigned char p[32], q[32],r[32],inv_r[32],inv_p[32],gf[256],fg[256],kkk[32]={0};

static const unsigned short t[256]={
  0,  1,128,255,154, 19,101,163, 53,173, 67, 62, 80, 93,170,115,197, 30,235, 87,214,243, 50, 65,149,  6, 79, 32, 70, 54,165, 24, 97,208,153,253,139,216, 56,143, 34,103, 63, 52,104,210,205,  8,244,241, 45,228,140, 29,171,199, 96, 64, 94, 41,135, 46,180,123,102,167, 15, 83,174, 35,126, 95,109,183, 99,191, 31,148, 91,231,176,225, 75,158,178,217, 69,155,138, 17, 20,144, 13,177,108,254,146, 51, 36,179,113,246, 74,202,192, 23,  2,120,198, 42,122, 10,230,185, 77, 61,212,219, 71, 86, 55,207,234,166,223,106,201,249,203, 39,156,  4,193,187, 14, 49,112,  7, 48,168,127,131, 84,118, 60,141, 68,152,125, 44, 40, 37, 25,226,116,121, 98,100,209, 16,233,107,252,161,186,196,  3,147,242,134,227,211,251,110,181, 85,237, 92, 88, 47,238,213, 66,142,218, 76,105,162,188,189,229,136, 89, 90,232,124,157,204,114,133,190,164,137,159, 58,111,160,200,215,195,117,239, 27, 43, 82,221,250, 38,220, 21, 73,172,119,182, 81,236,  5,132,151, 18,169, 33, 57,240,248,150,224, 12,184,194,222, 26, 22,247, 11,129, 28,130,206,175,145, 78, 59,  9, 72,245,};

static const unsigned short inv_t[256]={
  0,  1,106,166,131,226, 25,137, 47,253,111,244,237, 92,134, 66,159, 89,229,  5, 90,219,242,105, 31,152,241,212,246, 53, 17, 76, 27,231, 40, 69, 98,151,217,129,150, 59,109,213,149, 50, 61,179,138,135, 22, 97, 43,  8, 29,120, 38,232,204,252,144,115, 11, 42, 57, 23,182, 10,146, 86, 28,118,254,220,102, 82,185,114,251, 26, 12,224,214, 67,142,175,119, 19,178,192,193, 78,177, 13, 58, 71, 56, 32,156, 74,157,  6, 64, 41, 44,186,125,161, 94, 72,173,205,136,100,198, 15,154,210,143,222,107,155,110, 63,195,148, 70,140,  2,245,247,141,227,199,169, 60,191,202, 88, 36, 52,145,183, 39, 91,250, 96,167, 77, 24,235,228,147, 34,  4, 87,130,196, 83,203,206,163,187,  7,201, 30,123, 65,139,230, 14, 54,221,  9, 68,249, 80, 93, 84, 99, 62,174,223, 73,238,113,164,133,188,189,200, 75,104,132,239,209,165, 16,108, 55,207,126,103,128,197, 46,248,121, 33,158, 45,171,116,181, 20,208, 37, 85,184,117,218,215,240,124,236, 81,153,170, 51,190,112, 79,194,160,122, 18,225,176,180,211,233, 49,168, 21, 48,255,101,243,234,127,216,172,162, 35, 95,  3,};

/*
 * Number of columns (32-bit words) comprising the State. For this 
 * standard, Nb = 4.
 */
int Nb=4;

/*
 * Number of 32-bit words comprising the Cipher Key. For this 
 * standard, Nk = 4, 6, or 8.
 */
int Nk;

/*
 * Number of rounds, which is a function of  Nk  and  Nb (which is 
 * fixed). For this standard, Nr = 10, 12, or 14.
 */
int Nr;

/*
 * S-box transformation table
 */
static const uint8_t s_box[256] = {
	// 0     1     2     3     4     5     6     7     8     9     a     b     c     d     e     f
	0x63, 0x7c, 0x77, 0x7b, 0xf2, 0x6b, 0x6f, 0xc5, 0x30, 0x01, 0x67, 0x2b, 0xfe, 0xd7, 0xab, 0x76, // 0
	0xca, 0x82, 0xc9, 0x7d, 0xfa, 0x59, 0x47, 0xf0, 0xad, 0xd4, 0xa2, 0xaf, 0x9c, 0xa4, 0x72, 0xc0, // 1
	0xb7, 0xfd, 0x93, 0x26, 0x36, 0x3f, 0xf7, 0xcc, 0x34, 0xa5, 0xe5, 0xf1, 0x71, 0xd8, 0x31, 0x15, // 2
	0x04, 0xc7, 0x23, 0xc3, 0x18, 0x96, 0x05, 0x9a, 0x07, 0x12, 0x80, 0xe2, 0xeb, 0x27, 0xb2, 0x75, // 3
	0x09, 0x83, 0x2c, 0x1a, 0x1b, 0x6e, 0x5a, 0xa0, 0x52, 0x3b, 0xd6, 0xb3, 0x29, 0xe3, 0x2f, 0x84, // 4
	0x53, 0xd1, 0x00, 0xed, 0x20, 0xfc, 0xb1, 0x5b, 0x6a, 0xcb, 0xbe, 0x39, 0x4a, 0x4c, 0x58, 0xcf, // 5
	0xd0, 0xef, 0xaa, 0xfb, 0x43, 0x4d, 0x33, 0x85, 0x45, 0xf9, 0x02, 0x7f, 0x50, 0x3c, 0x9f, 0xa8, // 6
	0x51, 0xa3, 0x40, 0x8f, 0x92, 0x9d, 0x38, 0xf5, 0xbc, 0xb6, 0xda, 0x21, 0x10, 0xff, 0xf3, 0xd2, // 7
	0xcd, 0x0c, 0x13, 0xec, 0x5f, 0x97, 0x44, 0x17, 0xc4, 0xa7, 0x7e, 0x3d, 0x64, 0x5d, 0x19, 0x73, // 8
	0x60, 0x81, 0x4f, 0xdc, 0x22, 0x2a, 0x90, 0x88, 0x46, 0xee, 0xb8, 0x14, 0xde, 0x5e, 0x0b, 0xdb, // 9
	0xe0, 0x32, 0x3a, 0x0a, 0x49, 0x06, 0x24, 0x5c, 0xc2, 0xd3, 0xac, 0x62, 0x91, 0x95, 0xe4, 0x79, // a
	0xe7, 0xc8, 0x37, 0x6d, 0x8d, 0xd5, 0x4e, 0xa9, 0x6c, 0x56, 0xf4, 0xea, 0x65, 0x7a, 0xae, 0x08, // b
	0xba, 0x78, 0x25, 0x2e, 0x1c, 0xa6, 0xb4, 0xc6, 0xe8, 0xdd, 0x74, 0x1f, 0x4b, 0xbd, 0x8b, 0x8a, // c
	0x70, 0x3e, 0xb5, 0x66, 0x48, 0x03, 0xf6, 0x0e, 0x61, 0x35, 0x57, 0xb9, 0x86, 0xc1, 0x1d, 0x9e, // d
	0xe1, 0xf8, 0x98, 0x11, 0x69, 0xd9, 0x8e, 0x94, 0x9b, 0x1e, 0x87, 0xe9, 0xce, 0x55, 0x28, 0xdf, // e
	0x8c, 0xa1, 0x89, 0x0d, 0xbf, 0xe6, 0x42, 0x68, 0x41, 0x99, 0x2d, 0x0f, 0xb0, 0x54, 0xbb, 0x16};// f

/*
 * Inverse S-box transformation table
 */
static const uint8_t inv_s_box[256] = {
	// 0     1     2     3     4     5     6     7     8     9     a     b     c     d     e     f
	0x52, 0x09, 0x6a, 0xd5, 0x30, 0x36, 0xa5, 0x38, 0xbf, 0x40, 0xa3, 0x9e, 0x81, 0xf3, 0xd7, 0xfb, // 0
	0x7c, 0xe3, 0x39, 0x82, 0x9b, 0x2f, 0xff, 0x87, 0x34, 0x8e, 0x43, 0x44, 0xc4, 0xde, 0xe9, 0xcb, // 1
	0x54, 0x7b, 0x94, 0x32, 0xa6, 0xc2, 0x23, 0x3d, 0xee, 0x4c, 0x95, 0x0b, 0x42, 0xfa, 0xc3, 0x4e, // 2
	0x08, 0x2e, 0xa1, 0x66, 0x28, 0xd9, 0x24, 0xb2, 0x76, 0x5b, 0xa2, 0x49, 0x6d, 0x8b, 0xd1, 0x25, // 3
	0x72, 0xf8, 0xf6, 0x64, 0x86, 0x68, 0x98, 0x16, 0xd4, 0xa4, 0x5c, 0xcc, 0x5d, 0x65, 0xb6, 0x92, // 4
	0x6c, 0x70, 0x48, 0x50, 0xfd, 0xed, 0xb9, 0xda, 0x5e, 0x15, 0x46, 0x57, 0xa7, 0x8d, 0x9d, 0x84, // 5
	0x90, 0xd8, 0xab, 0x00, 0x8c, 0xbc, 0xd3, 0x0a, 0xf7, 0xe4, 0x58, 0x05, 0xb8, 0xb3, 0x45, 0x06, // 6
	0xd0, 0x2c, 0x1e, 0x8f, 0xca, 0x3f, 0x0f, 0x02, 0xc1, 0xaf, 0xbd, 0x03, 0x01, 0x13, 0x8a, 0x6b, // 7
	0x3a, 0x91, 0x11, 0x41, 0x4f, 0x67, 0xdc, 0xea, 0x97, 0xf2, 0xcf, 0xce, 0xf0, 0xb4, 0xe6, 0x73, // 8
	0x96, 0xac, 0x74, 0x22, 0xe7, 0xad, 0x35, 0x85, 0xe2, 0xf9, 0x37, 0xe8, 0x1c, 0x75, 0xdf, 0x6e, // 9
	0x47, 0xf1, 0x1a, 0x71, 0x1d, 0x29, 0xc5, 0x89, 0x6f, 0xb7, 0x62, 0x0e, 0xaa, 0x18, 0xbe, 0x1b, // a
	0xfc, 0x56, 0x3e, 0x4b, 0xc6, 0xd2, 0x79, 0x20, 0x9a, 0xdb, 0xc0, 0xfe, 0x78, 0xcd, 0x5a, 0xf4, // b
	0x1f, 0xdd, 0xa8, 0x33, 0x88, 0x07, 0xc7, 0x31, 0xb1, 0x12, 0x10, 0x59, 0x27, 0x80, 0xec, 0x5f, // c
	0x60, 0x51, 0x7f, 0xa9, 0x19, 0xb5, 0x4a, 0x0d, 0x2d, 0xe5, 0x7a, 0x9f, 0x93, 0xc9, 0x9c, 0xef, // d
	0xa0, 0xe0, 0x3b, 0x4d, 0xae, 0x2a, 0xf5, 0xb0, 0xc8, 0xeb, 0xbb, 0x3c, 0x83, 0x53, 0x99, 0x61, // e
	0x17, 0x2b, 0x04, 0x7e, 0xba, 0x77, 0xd6, 0x26, 0xe1, 0x69, 0x14, 0x63, 0x55, 0x21, 0x0c, 0x7d};// f



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

//#define Nb 8
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
 * Transformation in the Cipher that takes all of the columns of the 
 * State and mixes their data (independently of one another) to 
 * produce new columns.
 */
void mix_columns(uint8_t *state) {

	uint8_t a[] = {0x02, 0x01, 0x01, 0x03}; // a(x) = {02} + {01}x + {01}x2 + {03}x3
	uint8_t i, j, col[4], res[4];

	for (j = 0; j < Nb; j++) {
		for (i = 0; i < 4; i++) {
			col[i] = state[Nb*i+j];
		}

		coef_mult(a, col, res);

		for (i = 0; i < 4; i++) {
			state[Nb*i+j] = res[i];
		}
	}
}


/*
 * Transformation in the Inverse Cipher that is the inverse of 
 * MixColumns().
 */
void inv_mix_columns(uint8_t *state) {

	uint8_t a[] = {0x0e, 0x09, 0x0d, 0x0b}; // a(x) = {0e} + {09}x + {0d}x2 + {0b}x3
	uint8_t i, j, col[4], res[4];

	for (j = 0; j < Nb; j++) {
		for (i = 0; i < 4; i++) {
			col[i] = state[Nb*i+j];
		}

		coef_mult(a, col, res);

		for (i = 0; i < 4; i++) {
			state[Nb*i+j] = res[i];
		}
	}
}

/*
 * Transformation in the Cipher that processes the State using a non­
 * linear byte substitution table (S-box) that operates on each of the 
 * State bytes independently. 
 */
void sub_bytes(uint8_t *state) {

	uint8_t i, j;
	
	for (i = 0; i < 4; i++) {
		for (j = 0; j < Nb; j++) {
			// s_box row: yyyy ----
			// s_box col: ---- xxxx
			// s_box[16*(yyyy) + xxxx] == s_box[yyyyxxxx]
			state[Nb*i+j] = s_box[state[Nb*i+j]];
		}
	}
}

/*
 * Transformation in the Inverse Cipher that is the inverse of 
 * SubBytes().
 */
void inv_sub_bytes(uint8_t *state) {

	uint8_t i, j;

	for (i = 0; i < 4; i++) {
		for (j = 0; j < Nb; j++) {
			state[Nb*i+j] = inv_s_box[state[Nb*i+j]];
		}
	}
}

/*
 * Transformation in the Cipher that processes the State by cyclically 
 * shifting the last three rows of the State by different offsets. 
 */
void shift_rows(uint8_t *state) {

	uint8_t i, k, s, tmp;

	for (i = 1; i < 4; i++) {
		// shift(1,4)=1; shift(2,4)=2; shift(3,4)=3
		// shift(r, 4) = r;
		s = 0;
		while (s < i) {
			tmp = state[Nb*i+0];
			
			for (k = 1; k < Nb; k++) {
				state[Nb*i+k-1] = state[Nb*i+k];
			}

			state[Nb*i+Nb-1] = tmp;
			s++;
		}
	}
}

/*
 * Transformation in the Inverse Cipher that is the inverse of 
 * ShiftRows().
 */
void inv_shift_rows(uint8_t *state) {

	uint8_t i, k, s, tmp;

	for (i = 1; i < 4; i++) {
		s = 0;
		while (s < i) {
			tmp = state[Nb*i+Nb-1];
			
			for (k = Nb-1; k > 0; k--) {
				state[Nb*i+k] = state[Nb*i+k-1];
			}

			state[Nb*i+0] = tmp;
			s++;
		}
	}
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
    for(int i=0;i<32;i++)
    r[i]=q[i];
	for(int i=0;i<32;i++)
	inv_r[r[i]]=i;
}

void reverse(){
    for(int i=0;i<32;i++)
    q[i]=inv_p[r[p[i]]];
    for(int i=0;i<32;i++)
    r[i]=q[i];
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
        //n=c[i];
        //c[i]=inv_s_box[((n % 16) + (n >> 4) * 16)];
        c[i] = (256 + c[i] - k[r[i]]) % 256;
    }

}


void enc(uint8_t *m, __uint8_t *k)
{
    int i;
    unsigned char n=0,ff[32]={0};
	
    rounder();

		for(int l=0;l<32;l++)
		ff[l]=t[k[r[l]]];
		memcpy(k,ff,32);
    for (i = 0; i < 4; i++)
	{
		unsigned char tmp[32]={0};

		for(int j=0;j<8;j++){
		n= (m[(i*8+j)%32] + (ff[r[(i*8+j)%32]])) % 256;
		m[(i*8+j)]=t[((n % 16) + (n >> 4) * 16)];
		}

    }
    //shift_rows(m);
    //mix_columns(m);
}

void dec(uint8_t *c, uint8_t *k)
{
    int i;
    unsigned char n=0,ff[32]={0};

    //inv_mix_columns(c);
    //inv_shift_rows(c);
    for (i = 0; i < 4; i++){
		unsigned char tmp[32]={0};
		
		for(int j=0;j<8;j++){
		n=c[i*8+j];
	    c[i*8+j]=inv_t[((n % 16) + (n >> 4) * 16)];
		c[(i*8+j)%32] = (256 + c[(i*8+j)%32] - (ff[r[(i*8+j)%32]])) % 256;
		//
		}
    }

		for(int l=0;l<32;l++)
		ff[l]=inv_t[k[inv_r[l]]];
		memcpy(k,ff,32);

	    reverse();
}


/*
 * Performs the AES cipher operation
 */
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

/*
 * Performs the AES inverse cipher operation
 */
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

void t_box(unsigned *gf,unsigned *fg)
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
	//uint16_t t[256],inv_t[256];
	
	//for(i=0;i<256;i++){
	t_box(tt,inv_tt);
	//}
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

	aes_key_expansion(k, w);


    random_shuffle(p,SIZE_OF_ARRAY(p));
    random_shuffle(r,SIZE_OF_ARRAY(r));
	memcpy(s,r,32);
    for(i=0;i<32;i++){
        inv_p[p[i]]=i;
		inv_r[r[i]]=i;
    }

	printf("Plaintext message:\n");
	for (i = 0; i < 4; i++) {
		printf("%02x %02x %02x %02x ", m[4*i+0], m[4*i+1], m[4*i+2], m[4*i+3]);
	}

	printf("\n");
	for(i=0;i<32;i++)
	kkk[i]=255;
	//aes_cipher(m /* in */, out /* out */, w /* expanded key */);
	for(i=0;i<16;i++)
	enc(kkk,w);
	for(i=0;i<32;i++)
	m[i]^=kkk[i];
	printf("Ciphered message:\n");
	for (i = 0; i < 4; i++) {
		printf("%02x %02x %02x %02x ", m[4*i+0], m[4*i+1], m[4*i+2], m[4*i+3]);
		//printf("%02x %02x %02x %02x ", out[4*i+0], out[4*i+1], out[4*i+2], out[4*i+3]);
	}

	printf("\n");

	//aes_inv_cipher(out, m, w);
	//for(i=0;i<16;i++)
	//dec(kkk,w);
	memset(kkk,0,32);
	for(i=0;i<32;i++)
	kkk[i]=255;
	memcpy(r,s,32);
	memcpy(w,ss,32);
	for(i=0;i<16;i++)
	enc(kkk,w);
	
	for(i=0;i<32;i++)
	m[i]^=kkk[i];
	printf("\n");
	printf("Original message (after inv cipher):\n");
	for (i = 0; i < 4; i++) {
		printf("%02x %02x %02x %02x ", m[4*i+0], m[4*i+1], m[4*i+2], m[4*i+3]);
	}

	printf("\n");

	free(w);

}