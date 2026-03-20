#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <cblas.h>
#include <string.h>

#define dtype double
#define NITER 10
#define WARMUP 2

#include "include/my_time_lib.h"

int main(int argc, char *argv[]) {

    // Example $./gemm_test 1024 
    if (argc != 2) {
        printf("Usage: %s n \n\n"
               "Where:\n"
               "\tmatrix A has size n x n\n"
               "\tmatrix B has size n x n\n"
               "\tresult matrix C has size n x n\n",
               argv[0]);
        return 1;
    }

    dtype *A, *B, *C;
    double timers[NITER];

    /* Generate now two square matrices A and B of size (2^n) and fill them with random doubles.
     *  After this, compute the matrix multiplication C = A x B using OpenBLAS: cblas_dgemm()
     *
     * Tasks:
     *      1) Record CPU time on time[] vector. 
     *      2) Compute their arithmetic mean and geometric mean.
     *
     * NOTE:
     *      1) Fill A and B with random doubles
     *      2) Init C with 0.0
     */

    /* |========================================| */
    /* |           Put here your code           | */
    /* |========================================| */



    /*  Here we compute the vectors' arithmetic mean and geometric mean; 
     *  these functions must be implemented inside
     *   of the library "src/my_time_lib.c" (and their headers in "include/my_time_lib.h").
     */
    double a_mean = 0.0, g_mean = 0.0;


    /* |========================================| */
    /* |           Put here your code           | */
    /* |========================================| */


    printf(" %10s | %10s | %10s |\n", "v name", "arithmetic mean", "geometric mean");
    printf(" %10s | %10f | %10f |\n", "time", a_mean, g_mean);

    return(0);
}
