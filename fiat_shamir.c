#include <stdio.h>
#include "fiat_shamir.h"

// long long int p[10] = {268435455, 268435455, 268435455, 4095, 0, 0, 16777216, 0, 268435455, 15}; // 256-bit prime
// // p = 0xffffffff00000001000000000000000000000000ffffffffffffffffffffffff
// long long int g[10] = {144229014, 169055325, 187932916, 131601118, 15890179, 240532036, 133741798, 236110884, 186110450, 6};  // generator of Z*_p
// long long int p_1[10] = {268435454, 268435455, 268435455, 4095, 0, 0, 16777216, 0, 268435455, 15};  //p_1 = p-1




int main()
{
    /* Prover's part */
    long long int x[10] = { 223111552,
                            131469255,
                            18984424,
                            106798713,
                            158436973,
                            14876239,
                            212040987,
                            114476603,
                            134152580,
                            13134 };  // Hash(pass)
    long long int extended_x[20] = {0};
    for(int i = 0; i < 10; i++)
        extended_x[i] = x[i];
    long long int x_mod_p[10] = {0};
    // Barrett(x, x_mod_p, 20);
    Barrett(extended_x, x_mod_p, p, 20);

    long long int y[10] = {0};
    sqm_mod_p(g, x_mod_p, y, 10);   // y = g^x  (mod p)   is done
    for(int i = 0; i<10; i++)
        printf("\n y[%d]= %lld",i, y[i]);

/* y : P sends to V */
    //Prover choose a random number v 
    long long int v[10] = {20190170,
                           140767076,
                           152046517,
                           163557714,
                           89552916,
                           1848764,
                           78838121,
                           137173203,
                           3663413,
                           110 };
    long long int t[10] = {0};
    sqm_mod_p(g, v, t, 10);   // t = g^v  is done
/* t : P sends to V */


    /* Verifier has y and t. Now he will chhose a random number c  and sends to the prover */
    long long int c[10] = { 246941407,
                            92511354,
                            119248285,
                            186369831,
                            262066571,
                            191977025,
                            185213022,
                            100583089,
                            2104427,
                            110} ;
/* c : V sends to P */


    // Prover will calculate  r = v - c*x  mod(p-1)
    long long int product_cx[20] = {0}, cx[10];
    mult(c, x_mod_p, product_cx, 10);   // product_cx = c*x  is done
    Barrett(product_cx, cx, p_1, 20);   // here we need reduction in  mod (p-1)

    int flag = check(v, cx, 10);  // retruns  1 if v > cx,  2 if v < cx,  3 if v = cx
    long long int r[10] = {0};  // r = v - cx  mod(p-1)
    if((flag == 1) || (flag == 3))   // when v >= cx
        sub(v, cx, r, 10);
    if(flag == 2)   // when v < cx
    {
        long long int temp[10] = {0};
        sub(cx, v, temp, 10);   // temp = cx - v   is done
        sub(p_1, temp, r, 10);   // r = (p-1) - (cx-v) = v - cx + p-1  is done
    }
/* r : P sends to V */    
        

    /* Verifier's part */
    // compute : g^r*y^c  (mod p) and to authenticate match it with  t
    long long int g_r[10] = {0}, y_c[10] = {0}, pr_gr_yc[20] = {0}, gr_yc[10]={0} ;
    sqm_mod_p(g, r, g_r, 10);   // g_r = g^r  (mod p)   is computed
    sqm_mod_p(y, c, y_c, 10);   // y_c = y^c  (mod p)   is computed
    mult(g_r, y_c, pr_gr_yc, 10);   // pr_gr_yc = g^r * y^c  is computed
    Barrett(pr_gr_yc, gr_yc, p, 20);  // gr_yc  =  pr_gr_yc (mod p)  is computed


/* Authentication part */
    int AUT = 0;
    for(int i = 0; i < 10; i++ )
    {
        if(t[i] != gr_yc[i])
        {
            printf("\n FAILED ! \n");
            break;
        }
        else
            AUT++;
        // printf("\n t[%d] = %lld              gr_yc[%d] = %lld  ",i, t[i],i, gr_yc[i]);
    }
    if(AUT == 10)
    printf("\n Authentication successful \n");

    return 0;
}
