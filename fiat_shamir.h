#ifndef FIAT_SHAMIR_DOT_H
#define FIAT_SHAMIR_DOT_H

extern long long int p[10];
extern long long int g[10];
extern long long int p_1[10];

void add(long long int a[], long long int b[], long long int c[], int size);  // c = a + b
void mult(long long int a[], long long int b[], long long int m[], int size); // m = a * b
void sub(long long int a[], long long int b[], long long int s[], int size);  // s = a - b  with assumption a > b
int check(long long int x[], long long int y[], int size); // returns 1 when x>y, 2 when x<y , 3 when x=y

void Barrett(long long int x[], long long int x_Barrett_p[], long long int n[], int size);    // Barrett reduction with ( mod p), where T = mu(p) = floor(B^2n/p),  B = 2^28, n = block size = 10
void sqm_mod_p(long long int y[], long long int d[], long long int z[], int size);    // z = y^b (mod p), size = 10,here


#endif