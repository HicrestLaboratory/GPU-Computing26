#include "../include/my_time_lib.h"

// Put here the implementation of arithmetic_mean and geometric_mean

double arithmetic_mean(double *v, int len) {

    double mu = 0.0;
    for (int i=0; i<len; i++)
        mu += (double)v[i];
    mu /= (double)len;

    return(mu);
}

double geometric_mean(double *v, int len) {

    double mu = 1.0;
    for (int i=0; i<len; i++) {
        mu *= (v[i] > 0) ? ((double)v[i]) : 1;
    }
    mu = pow(mu, 1.0 / len);

    return(mu);
}
      

// -------------------------------------------------
