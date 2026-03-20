#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <string.h>

#define dtype int
#define NITER 10
#define WARMUP 2

#include "include/my_time_lib.h"

void prefix_sum(dtype *A, dtype *P, size_t n)
{
    P[0] = A[0];

    for (size_t i = 1; i < n; i++)
        P[i] = P[i-1] + A[i];
}


int main(int argc, char *argv[]) {

    // Example $./prefixsum 32 
    if (argc != 2) {
        printf("Usage: %s n \n\n"
               "Where:\n"
               "\tarray A has size n\n"
               "\tresult array P has size n\n",
               argv[0]);
        return 1;
    }

    dtype *A, *P;
    double timers[NITER];

    /* Ex1: Generate now one array A and fill it with random ints.
     *  After this, compute the prefix sum P 
     *
     * Tasks:
     *      1) Record CPU time on time[] vector. 
     *      2) Compute their arithmetic mean and geometric mean.
     *
     * NOTE:
     *      1) Fill A with random int
     *      2) You could print input and result for small arrays to check the correctness
     */

    int n = atoi(argv[1]);
    A = (dtype*)malloc(n*sizeof(dtype));
    P = (dtype*)malloc(n*sizeof(dtype));

    if (!A || !P) {
        printf("Memory allocation failed\n");
        return 1;
    }

    /* initialize array */
    srand(0);
    for (size_t i=0;i<n;i++) {
        A[i] = rand()%10;
    }
        

    double iter_time;

    TIMER_DEF(0);

    for (int i=-WARMUP; i<NITER; i++) {

        TIMER_START(0);
        prefix_sum(A, P, n);
        TIMER_STOP(0);

        iter_time = TIMER_ELAPSED(0) / 1.e6;
        if( i >= 0) timers[i] = iter_time;
    }


    /* Ex1:  Here we compute the vectors' arithmetic mean and geometric mean; 
     *  these functions must be implemented inside
     *   of the library "src/my_time_lib.c" (and their headers in "include/my_time_lib.h").
     */
    double a_mean = 0.0, g_mean = 0.0;
    a_mean = arithmetic_mean(timers, NITER);
    g_mean = geometric_mean(timers, NITER);


    printf(" %10s | %10s | %10s |\n", "v name", "arithmetic mean", "geometric mean");
    printf(" %10s | %10f | %10f |\n", "time", a_mean, g_mean);

    /* Ex2: By varying the input array size, and data type, e.g., float, double, 
     *      Benchmark the effective bandwidth
     * Note:Bandwidth unit GB/s
     */
    double bytes = 0.0;
    double bandwidth = 0.0; 

    bytes = 2.0 * n * sizeof(dtype);
    bandwidth = bytes / a_mean / 1e9;

    printf(" %10s | %10s | %10s | %10s |\n", "array size","bytes","time(s)","bandwidth(GB/s)");
    printf(" %10s | %10f | %10f | %10f |\n", argv[1], bytes, a_mean, bandwidth);


    return(0);
}
