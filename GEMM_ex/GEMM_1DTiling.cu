#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <time.h>
#include <sys/time.h>
#include <math.h>
#include "include/helper_cuda.h"
#include <cuda_runtime.h>
#include <cub/cub.cuh>

#define TIMER_DEF     struct timeval temp_1, temp_2

#define TIMER_START   gettimeofday(&temp_1, (struct timezone*)0)

#define TIMER_STOP    gettimeofday(&temp_2, (struct timezone*)0)

#define TIMER_ELAPSED ((temp_2.tv_sec-temp_1.tv_sec)+(temp_2.tv_usec-temp_1.tv_usec)/1000000.0)

#define DBG_CHECK { printf("DBG_CHECK: file %s at line %d\n", __FILE__, __LINE__ ); }
#define DEBUG  // without debug (with random imputs) the kernel does not work

#define NPROBS 3
#define STR(s) #s
#define XSTR(s) STR(s)
#define dtype float

#define RUN_SOLUTIONS

#define PRINT_MATRIX(A, N, M, ST ) {  \
      int i, j;  \
      printf("%s:\n", ( ST ));  \
      for (i=0; i< ( N ); i++) {  \
        printf("\t");  \
        for (j=0; j< ( M ); j++)  \
          printf("%6.3f ", A[i*( M ) + j]);  \
        printf("\n");  \
      }  \
      printf("\n\n");  \
}

float matrix_error (int n, int m, const dtype* A, const dtype* B) {
  int i, j;
  dtype error = (dtype)0;
  for (i=0; i<n; i++)
    for (j=0; j<m; j++)
      error += fabs(B[i*m + j] - A[i*m + j]);

  return(error);
}

#define CEIL_DIV( N, D ) ((( N ) % ( D )) == 0) ? (( N )/( D )) : ((( N )/( D ))+1)



// --------------------- sgemm_block_tiling

#define BN 64
#define BK 16
#define BM 64

#define TM 4



// ----------------------- sgemm_warptiling

#define BNWARP 128
#define BKWARP 8
#define BMWARP 128
#define TMWARPS 2
#define TNWARPS 2

#define WARPNUM 32
#define WARPSN 4
#define WARPSM 8
#define WN 32
#define WM 16

#define WARPSIZE 32
#define WSUBN 8
#define WSUBM 4
#define WNITER 2
#define WMITER 2

// ----------------------------------------


int verbose;


__global__ void sgemm_block_tiling(int N, int K, int M, float alpha, const dtype *A, const dtype *B, float beta, dtype *C) {

  /* |========================================| */
  /* |           Put here your code           | */
  /* |========================================| */



}





__global__ void sgemm_warptiling(int N, int K, int M, float alpha, const dtype *A, const dtype *B, float beta, dtype *C) {

  /* |========================================| */
  /* |           Put here your code           | */
  /* |========================================| */
}

dtype* execute_gemm_kernel (int n, int k, int m, float alpha, dtype* A, dtype* B, float beta, void (*gemm_kernel)(int, int, int, float, const dtype*, const dtype*, float, dtype*), float* Bandwidth, float* CompTime, double* Flops) {
  int grd_sizeX, grd_sizeY;
  int blk_sizeX, blk_sizeY;

  // ---------------------------------
 if (gemm_kernel == sgemm_block_tiling) {
    if ( ((BN * BM) / TM) != (BN * BK) ) {
      fprintf(stderr, "[%d] ERROR: ((BN * BM) / TM) = %d != %d = (BN * BK)\n", __LINE__, ((BN * BM) / TM), (BN * BK));
      exit(42);
    }

    grd_sizeX = CEIL_DIV(n, BN);
    grd_sizeY = CEIL_DIV(m, BM);

    blk_sizeX = (BN * BM) / TM;
    blk_sizeY = 1;
  } else {
      if (gemm_kernel == sgemm_warptiling) {
        grd_sizeX = CEIL_DIV(n, BNWARP);
        grd_sizeY = CEIL_DIV(m, BMWARP);

        blk_sizeX = BNWARP * BKWARP;
        blk_sizeY = 1;

        unsigned int flag = 0U, e = 0U;
        (WARPNUM * 32 == blk_sizeX)           ? (e &= (1U<<0)) : (flag &= (1U<<0));

        ((WARPSN * WARPSM) == WARPNUM)        ? (e &= (1U<<1)) : (flag &= (1U<<1));
        ((WSUBN * WNITER) == WN)              ? (e &= (1U<<2)) : (flag &= (1U<<2));
        ((WN * WARPSN) == BNWARP)             ? (e &= (1U<<3)) : (flag &= (1U<<3));

        ((WSUBM * WMITER) == WM)              ? (e &= (1U<<4)) : (flag &= (1U<<4));
        ((WM * WARPSM) == BMWARP)             ? (e &= (1U<<5)) : (flag &= (1U<<5));

        if (flag != 0U) {
          fprintf(stderr, "ERROR: line %d, e = %u, flag = %u\n", __LINE__, e, flag);
          fprintf(stderr, "ERROR: WARPSN = %d, WSUBN = %d, WNITER = %d, WN = %d, WARPSM = %d, WSUBM = %d, WMITER = %d, WM = %d\n", WARPSN, WSUBN, WNITER, WN, WARPSM, WSUBM, WMITER, WM);
          exit(42);
        }
      }
    }
  }
  // ---------------------------------

  // ------------------- allocating GPU vectors ----------------------
  dtype *dev_A, *dev_B, *dev_C;

  checkCudaErrors( cudaMalloc(&dev_A, n*k*sizeof(dtype)) );
  checkCudaErrors( cudaMalloc(&dev_B, k*m*sizeof(dtype)) );
  checkCudaErrors( cudaMalloc(&dev_C, n*m*sizeof(dtype)) );
  size_t bandwidth_numerator = ((n*k) + (k*m) + (n*m))*sizeof(dtype);

  // ----------------- copy date from host to device -----------------

  checkCudaErrors( cudaMemcpy(dev_A, A, n*k*sizeof(dtype), cudaMemcpyHostToDevice) );
  checkCudaErrors( cudaMemcpy(dev_B, B, k*m*sizeof(dtype), cudaMemcpyHostToDevice) );
  checkCudaErrors( cudaMemset(dev_C, 0, n*m*sizeof(dtype)) );

  // ---------- compute GPU_tmp_b with the reduction kernel ----------
  TIMER_DEF;
  TIMER_START;

  {
      dim3 block_size(blk_sizeX, blk_sizeY, 1);
      dim3 grid_size(grd_sizeX, grd_sizeY, 1);
      printf("%d: block_size = (%d, %d), grid_size = (%d, %d)\n", __LINE__, block_size.x, block_size.y, grid_size.x, grid_size.y);
      gemm_kernel<<<grid_size, block_size>>>(n, k, m, alpha, (const dtype*)dev_A, (const dtype*)dev_B, beta, dev_C);
  }


  checkCudaErrors( cudaDeviceSynchronize() );
  TIMER_STOP;
  *CompTime += TIMER_ELAPSED;
  *Bandwidth = bandwidth_numerator / ((*CompTime)*1e+9);
  *Flops  = (n) / ((*CompTime)*1e+9);
  *Flops *= m;
  *Flops *= k;
  *Flops *= 2;

  // --------------- copy results from device to host ----------------

  dtype *GPU_C = (dtype*)malloc(sizeof(dtype)*n*m);
  checkCudaErrors( cudaMemcpy(GPU_C, dev_C, n*m*sizeof(dtype), cudaMemcpyDeviceToHost) );

  if (verbose > 1)
    PRINT_MATRIX(GPU_C, n, m, "GPU_C form execute_gemm_kernel")

  checkCudaErrors( cudaFree(dev_A) );
  checkCudaErrors( cudaFree(dev_B) );
  checkCudaErrors( cudaFree(dev_C) );

  return(GPU_C);
}


int main(int argc, char *argv[]) {

  printf("====================================== Problem computations ======================================\n");
  // =========================================== Set-up the problem ============================================

  if (argc < 3) {
    printf("Usage: lab3_ex1 e v [CPU_ON = 1]\n");
    return(1);
  }
  printf("argv[1] = %s\n", argv[1]);
  printf("argv[2] = %s\n", argv[2]);
  if (argc > 3)
    printf("argv[3] = %s\n", argv[3]);

  // ---------------- set-up the problem size -------------------


  int e = atoi(argv[1]), n = (1<<(e/2)), k = n, m = n, i, j, CPU_ON = 1;
  float alpha = 1.0f, beta = 1.0f;
  verbose = atoi(argv[2]);
  if (argc > 3)
    CPU_ON = atoi(argv[3]);

  // BUG check: code to generalize
  if ((e%2) != 0) {
    printf("Now the code only support squared matrices. So, since the generated matrix will have dimensions 2^(e/2) x 2^(e/2), e must be even\n");
    exit(42);
  }

  printf("e = %d --> n = k = m = 2^(e/2) = %d\n", e, n);
  printf("alpha = %f, beta = %f\n", alpha, beta);
  printf("CPU_ON = %d\n", CPU_ON);
  printf("verbose = %d\n", verbose);
  printf("dtype = %s\n", XSTR(dtype));

  // ======================================== Get the device properties ========================================
  printf("======================================= Device properties ========================================\n");

  int deviceCount = 0;
  cudaError_t error_id = cudaGetDeviceCount(&deviceCount);

  int dev;
  for (dev = 0; dev < deviceCount; ++dev) {
    cudaSetDevice(dev);
    cudaDeviceProp deviceProp;
    cudaGetDeviceProperties(&deviceProp, dev);

    printf("\nDevice %d: \"%s\"\n", dev, deviceProp.name);

    printf("  Memory Clock rate:                             %.0f Mhz\n",
           deviceProp.memoryClockRate * 1e-3f);

    printf("  Memory Bus Width:                              %d bit\n",
           deviceProp.memoryBusWidth);

    printf("  Peak Memory Bandwidth:                     %7.3f GB/s\n",
           2.0*deviceProp.memoryClockRate*(deviceProp.memoryBusWidth/8)/1.0e6);

    printf("  (%03d) Multiprocessors, (%03d) CUDA Cores/MP:    %d CUDA Cores\n",
           deviceProp.multiProcessorCount,
           _ConvertSMVer2Cores(deviceProp.major, deviceProp.minor),
           _ConvertSMVer2Cores(deviceProp.major, deviceProp.minor) *
               deviceProp.multiProcessorCount);

    printf("  Peak Arithmetic Intensity:                     %7.3f GFLOPS/s\n",
           2.0*deviceProp.memoryClockRate*(_ConvertSMVer2Cores(deviceProp.major, deviceProp.minor) *
               deviceProp.multiProcessorCount)/1.0e6);
    printf("  Shared Memory per SM:                     %d B\n", deviceProp.sharedMemPerMultiprocessor);
  }

  // ------------------ set-up the timers ---------------------

  TIMER_DEF;
  const char* lables[NPROBS] = {"CPU check", "Kernel 1", "Kernel 2"};
  float errors[NPROBS], Times[NPROBS], Bandwidths[NPROBS], error;
  double Flops[NPROBS];
  for (i=0; i<NPROBS; i++) {
    errors[i] = 1<<30;
    Bandwidths[i] = 0;
    Flops[i] = 0;
    Times[i] = 0;
  }


  // ------------------- set-up the problem -------------------

  dtype *A, *B, *GPU_C, *CPU_C;
  A = (dtype*)malloc(sizeof(dtype)*n*k);
  B = (dtype*)malloc(sizeof(dtype)*k*m);
  CPU_C = (dtype*)malloc(sizeof(dtype)*n*m);
//   GPU_C = (dtype*)malloc(sizeof(dtype)*n*m);

  time_t t;
  srand((unsigned) time(&t));


  for (i=0; i<(n*k); i++)
    A[i] = ((dtype)(i/m)/(dtype)m) + 1.0f;
  for (i=0; i<(k*m); i++)
    B[i] = (dtype)(1);

#ifdef DEBUG
  if (verbose > 0) {
    PRINT_MATRIX(A, n, k, "A")

    PRINT_MATRIX(B, k, m, "B")
  }
#endif
  // ======================================== Running the computations =========================================

  /* [ ... ]
   */

  // ========================== CPU computation =========================
  if (CPU_ON) {

    TIMER_START;
    for (i=0; i<n; i++)
      for (j=0; j<m; j++)
        for (int h=0; h<k; h++)
          CPU_C[i*m +j] += A[i*k + h] * B[h*m + j];
    TIMER_STOP;

    Times[0] = TIMER_ELAPSED;
    errors[0] = 0.0f;
    Bandwidths[0] = 0.0f;
    Flops[0]  = (n) / (Times[0]*1e+9);
    Flops[0] *= m;
    Flops[0] *= k;


    if (verbose > 0)
      PRINT_MATRIX(CPU_C, n, m, "CPU_C")

  } else {
    Times[0] = -1.0f;
    errors[0] = -1.0f;
    Bandwidths[0] = -1.0f;
    Flops[0] = -1.0f;
  }
  
  printf("=========================== GPU Kernel 1 ===========================\n");
  // =========================== GPU Kernel 1 ===========================

  GPU_C = execute_gemm_kernel(n, k, m, alpha, A, B, beta, sgemm_block_tiling, &Bandwidths[1], &Times[1], &Flops[1]);

  // ------------- Compare GPU and CPU solution --------------

  (CPU_ON) ? (error = matrix_error(n, m, CPU_C, GPU_C)) : (error = 0.0f) ;
  errors[1] = error;

  if (verbose > 0)
    PRINT_MATRIX(GPU_C, n, m, "GPU_C")

  free(GPU_C);



  printf("=========================== GPU Kernel 2 ===========================\n");
  // =========================== GPU Kernel 2 ===========================

  GPU_C = execute_gemm_kernel(n, k, m, alpha, A, B, beta, sgemm_warptiling, &Bandwidths[2], &Times[2], &Flops[2]);

  // ------------- Compare GPU and CPU solution --------------

  (CPU_ON) ? (error = matrix_error(n, m, CPU_C, GPU_C)) : (error = 0.0f) ;
  errors[2] = error;

  if (verbose > 0)
    PRINT_MATRIX(GPU_C, n, m, "GPU_C")

  free(GPU_C);


  printf("\n\n");
  if (!(CPU_ON)) printf("CPU check not lunched!!\n");
  printf("Solution\n %9s\t%9s\t%9s\t%16s\t%16s\n", "type", "error", "time (s)", "flops (GFLOPS/s)", "bandwidth (GB/s)");
  for (int i=0; i<NPROBS; i++) {
    if ((i != 6))
      printf("%12s:\t%9.6f\t%9.6f\t%16.6lf\t%16.6f\n", lables[i], errors[i], Times[i], Flops[i], Bandwidths[i]);
  }
  printf("\n");

  printf("GPU times: e Kernel1_time Kernel1_flops Kernel2_time Kernel2_flops ... on stderr\n");
  fprintf(stderr, "%d\t", e);
  for (i=1; i<NPROBS; i++)
    fprintf(stderr, "%f\t%f\t", Times[i], Flops[i]);
  fprintf(stderr, "\n");

  return(0);
}
