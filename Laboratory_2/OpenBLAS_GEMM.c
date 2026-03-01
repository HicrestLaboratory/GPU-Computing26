#include <stdio.h>
#include <cblas.h>

int main() {
    int m = 2, n = 2, k = 2;
    double A[4] = {1.0, 2.0, 3.0, 4.0};
    double B[4] = {5.0, 6.0, 7.0, 8.0};
    double C[4] = {0.0, 0.0, 0.0, 0.0};

    // Perform C = A * B using BLAS
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                m, n, k, 1.0, A, k, B, n, 0.0, C, n);

    // Print result
    printf("Result matrix C:\n");
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            printf("%lf ", C[i * n + j]);
        }
        printf("\n");
    }
    return 0;
}
