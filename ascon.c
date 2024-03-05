#include <stdio.h>

long long state[5] = { 0 }, t[5] = { 0 };
long long constants[16] = {0xf0, 0xe1, 0xd2, 0xc3, 0xb4, 0xa5, 0x96, 0x87, 0x78, 0x69, 0x5a, 0x4b, 0x3c, 0x2d, 0x1e, 0x0f};

long long print_state(long long state[5]){
   for(int i = 0; i < 5; i++){
      printf("%llx\n", state[i]);
   } 
}

long long rotate(long long x, int l) {
   long long temp;
   temp = (x >> l) ^ (x << (64 - l));
   return temp;
}

void add_constant(long long state[5], int i, int a) {
   // Menambah konstan pada state blok ke 2 sesuai dengan spec Ascon
   state[2] = state[2] ^ constants[12 - a + i];
}
void sbox(long long x[5]) {
   // Mensubtitusikan angka menjadi angka baru pada state sesuai dengan sbox
   // Biasanya sbox dilakukan dengan menggunakan tabel lookup
   // tapi juga bisa menggunakan bitslice implementation sesuai dengan spec ascon
   // justru bitwise operation lebih bagus karena lebih ringan
   // dan juga menghindari penggunaan side channel attack.

   x[0] ^= x[4]; x[4] ^= x[3]; x[2] ^= x[1];
   t[0] = x[0]; t[1] = x[1]; t[2] = x[2]; t[3] = x[3]; t[4] = x[4];
   t[0] =~ t[0]; t[1] =~ t[1]; t[2] =~ t[2]; t[3] =~ t[3]; t[4] =~ t[4];
   t[0] &= x[1]; t[1] &= x[2]; t[2] &= x[3]; t[3] &= x[4]; t[4] &= x[0];
   x[0] ^= t[1]; x[1] ^= t[2]; x[2] ^= t[3]; x[3] ^= t[4]; x[4] ^= t[0];
   x[1] ^= x[0]; x[0] ^= x[4]; x[3] ^= x[2]; x[2] =~ x[2];
}
void linear(long long state[5]) {
   // Kita akan melakukan operasi rotasi terhadap state dengan tiap
   // 64 bit memiliki rotasi yang berbeda.
   // besar bit rotasi ditentukan pada spec ascon paper.
   
   long long temp0, temp1;
   temp0 = rotate(state[0], 19);
   temp1 = rotate(state[0], 28);
   state[0] ^= temp0 ^ temp1;
   temp0 = rotate(state[1], 61);
   temp1 = rotate(state[1], 39);
   state[1] ^= temp0 ^ temp1;
   temp0 = rotate(state[2], 1);
   temp1 = rotate(state[2], 6);
   state[2] ^= temp0 ^ temp1;
   temp0 = rotate(state[3], 10);
   temp1 = rotate(state[3], 17);
   state[3] ^= temp0 ^ temp1;
   temp0 = rotate(state[4], 7);
   temp1 = rotate(state[4], 41);
   state[4] ^= temp0 ^ temp1;
}

void p(long long state[5], int a){
   for (int i = 0; i < a; i++){
      add_constant(state, i, a);
      sbox(state);
      linear(state);
   }
}

void initialization(long long state[5], long long key[2]) {
   p(state, 12);
   state[3] ^= key[0];
   state[4] ^= key[1];
}

void associated_data(long long state[5], int length, long long associated_data_text[]) {
   for (int i = 0; i < length; i++){
      state[0] = associated_data_text[i] ^ state[0];
      p(state, 6);
   }
   state[5] = state[5] ^ 0x0000000000000001;
}

void finalization(long long state[5], long long key[2]) {
   state[1] ^= key[0];
   state[2] ^= key[1];
   p(state, 12);
   state[3] ^= key[0];
   state[4] ^= key[1];

}

void encrypt(long long state[5], int length, long long plaintext[], long long ciphertext[]) {
   ciphertext[0] = plaintext[0] ^ state[0];
   for (int i = 1; i < length; i++){
      p(state, 6);
      ciphertext[i] = plaintext[i] ^ state[0];
      state[0] = ciphertext[i];
   }
}

void decrypt(long long state[5], int length, long long plaintext[], long long ciphertext[]){
   plaintext[0] = ciphertext[0] ^ state[0];
   for (int i = 1; i < length; i++){
      p(state, 6);
      plaintext[i] = ciphertext[i] ^ state[0];
      state[0] = ciphertext[i];
   }
}


int main() {
   // initialize nonce, key and IV
   long long nonce[2] = { 0x0000000000000001, 0x0000000000000002 };
   long long key[2] = { 0 };
   long long IV = 0x80400c0600000000;
   long long plaintext[] = {0x1234567890abcdef, 0x1234567890abcdef};
   long long ciphertext[2] = { 0 };
   long long associated_data_text[] = { 0x787878, 0x878787, 0x09090};

   //encryption
   //initialize state
   state[0] = IV;
   state[1] = key[0];
   state[2] = key[1];
   state[3] = nonce[0];
   state[4] = nonce[1];
   initialization(state,key);
   associated_data(state, 3, associated_data_text);
   print_state(state);
   encrypt(state, 2, plaintext, ciphertext);
   printf("\nciphertext:%llx%llx\n", ciphertext[0], ciphertext[1]);
   finalization(state, key);
   printf("tag:%llx%llx\n", state[3], state[4]);



   //decryption
        
   long long ciphertextdecrypt[2] = {0};
   for(int i = 0; i < 2; i++){
      ciphertextdecrypt[i] = ciphertext[i];
   }
   long long plaintextdecrypt[10] = { 0 };

   //initialize state
   state[0] = IV;
   state[1] = key[0];
   state[2] = key[1];
   state[3] = nonce[0];
   state[4] = nonce[1];

   initialization(state,key);
   print_state(state);
   associated_data(state, 3, associated_data_text);
   decrypt(state, 2, plaintextdecrypt, ciphertextdecrypt);
   printf("\nplaintext:%llx%llx\n", plaintextdecrypt[0], plaintextdecrypt[1]);
   finalization(state, key);
   printf("tag:%llx%llx\n", state[3], state[4]);}