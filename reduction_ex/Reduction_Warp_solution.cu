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

#define NPROBS 2
#define ITER 10
#define STR(s) #s
#define XSTR(s) STR(s)
#define dtype float

#define BLK_SIZE 256

int verbose;


__inline__ __device__ dtype warpReduceSum(dtype val) {
    for (int offset = warpSize/2; offset > 0; offset /= 2)
        val += __shfl_down_sync(0xffffffff, val, offset);
    return val;
}

__inline__ __device__ dtype blockReduceSum(dtype val) {
    __shared__ dtype shared[32]; // max 1024 threads / 32 = 32 warps

    int lane = threadIdx.x % warpSize;
    int wid  = threadIdx.x / warpSize;

    val = warpReduceSum(val);  // each warp performs partial reduction

    if (lane == 0)
        shared[wid] = val;    // write reduced value to shared memory

    __syncthreads();

    val = (threadIdx.x < blockDim.x / warpSize) ? shared[lane] : 0;

    if (wid == 0)
        val = warpReduceSum(val); //Final reduce within first warp

    return val;
}

__global__ void kernel_reduction_warp(int len, dtype* g_idata, dtype* g_odata) {

    unsigned int tid = blockIdx.x * blockDim.x + threadIdx.x;

    dtype sum = 0;

    // grid-stride loop: reduce multiple elements per thread
    for (int i = tid; i < len; i += blockDim.x * gridDim.x) {
        sum += g_idata[i];
    }

    sum = blockReduceSum(sum);

    // one value per block 
    if (threadIdx.x == 0) {
        g_odata[blockIdx.x] = sum;
    }
}



dtype execute_reduction_kernel (int len, dtype* a, void(*reduction_kernel)(int,dtype*,dtype*), float* Bandwidth, float* CompTime) {
  int grd_size = len/BLK_SIZE;
  int kernel_len = len;
  dtype *GPU_tmp_b;

  // ------------------- allocating GPU vectors ----------------------
  dtype *dev_a, *dev_tmp_b, *dev_b;

  checkCudaErrors( cudaMalloc(&dev_a, len*sizeof(dtype)) );
  checkCudaErrors( cudaMalloc(&dev_tmp_b, grd_size*sizeof(dtype)) );
  checkCudaErrors( cudaMalloc(&dev_b, sizeof(dtype)) );
  size_t bandwidth_numerator = 0;

  // ----------------- copy date from host to device -----------------

  checkCudaErrors( cudaMemcpy(dev_a, a, len*sizeof(dtype), cudaMemcpyHostToDevice) );
  checkCudaErrors( cudaMemset(dev_tmp_b, 0, grd_size*sizeof(dtype)) );
  checkCudaErrors( cudaMemset(dev_b, 0, sizeof(dtype)) );

  // -------- Warm-up runs --------
  for (int w = 0; w < ITER; w++) {
      dim3 block_size(BLK_SIZE, 1, 1);
      dim3 grid_size(grd_size, 1, 1);

      reduction_kernel<<<grid_size, block_size>>>(kernel_len, dev_a, dev_tmp_b);
  }
  cudaDeviceSynchronize();


  // ---------- Fisr reduction pass: compute GPU_tmp_b with the reduction kernel ----------
  TIMER_DEF;
  TIMER_START;

  {

      dim3 block_size(BLK_SIZE, 1, 1);
      dim3 grid_size(grd_size, 1, 1);
      printf("%d: block_size = %d, grid_size = %d, kernel_len = %d\n", __LINE__, block_size.x, grid_size.x, kernel_len);
      reduction_kernel<<<grid_size, block_size>>>(kernel_len, dev_a, dev_tmp_b);
      bandwidth_numerator += kernel_len*sizeof(dtype);
      kernel_len = grd_size;
      bandwidth_numerator += kernel_len*sizeof(dtype);
  }


  checkCudaErrors( cudaDeviceSynchronize() );
  TIMER_STOP;
  *CompTime += TIMER_ELAPSED;

  // ------------------------ debug if needed --------------------
  if (verbose > 0) {
    GPU_tmp_b = (dtype*) malloc(sizeof(dtype)*kernel_len);
    checkCudaErrors( cudaMemcpy(GPU_tmp_b, dev_tmp_b, kernel_len*sizeof(dtype), cudaMemcpyDeviceToHost) );
    dtype CPUsum_of_GPU_tmp_b = 0.0f;

    if (verbose > 1) printf("GPU_tmp_b: ");
    for (int i=0; i<kernel_len; i++) {
        if (verbose > 1) {
          if (i==0 || GPU_tmp_b[i] != GPU_tmp_b[i-1])
            printf("%f ... ", GPU_tmp_b[i]);
        }
        CPUsum_of_GPU_tmp_b += GPU_tmp_b[i];
    }
    if (verbose > 1) printf("\n");
    printf("CPUsum_of_GPU_tmp_b = %f\n", CPUsum_of_GPU_tmp_b);

    free(GPU_tmp_b);
  }

  // --------------------- Multi-pass reduction ---------------------------------
  while (grd_size > BLK_SIZE) {
    dtype* dev_tmp_b2;
    checkCudaErrors( cudaMalloc(&dev_tmp_b2, sizeof(dtype)*(kernel_len/BLK_SIZE)) );
    checkCudaErrors( cudaMemset(dev_tmp_b2, 0, sizeof(dtype)*(kernel_len/BLK_SIZE)) );
    // ---------- compute GPU_tmp_b with the reduction kernel ----------
    TIMER_START;

    {
        dim3 block_size(BLK_SIZE, 1, 1);
        grd_size = kernel_len / BLK_SIZE;
        dim3 grid_size(grd_size, 1, 1);
        printf("%d: block_size = %d, grid_size = %d, kernel_len = %d\n", __LINE__, block_size.x, grid_size.x, kernel_len);
        reduction_kernel<<<grid_size, block_size>>>(kernel_len, dev_tmp_b, dev_tmp_b2);
        bandwidth_numerator += kernel_len*sizeof(dtype);
        kernel_len = grd_size;
        bandwidth_numerator += kernel_len*sizeof(dtype);
    }


    checkCudaErrors( cudaDeviceSynchronize() );
    TIMER_STOP;
    *CompTime += TIMER_ELAPSED;
    checkCudaErrors( cudaFree(dev_tmp_b) );
    dev_tmp_b = dev_tmp_b2;

    // ------------------------ debug if needed --------------------
    if (verbose > 0) {
        GPU_tmp_b = (dtype*) malloc(sizeof(dtype)*kernel_len);
        checkCudaErrors( cudaMemcpy(GPU_tmp_b, dev_tmp_b, kernel_len*sizeof(dtype), cudaMemcpyDeviceToHost) );
        dtype CPUsum_of_GPU_tmp_b = 0.0f;

        if (verbose > 1) printf("GPU_tmp_b: ");
        for (int i=0; i<kernel_len; i++) {
          if (verbose > 1) {
            if (i==0 || GPU_tmp_b[i] != GPU_tmp_b[i-1])
              printf("%f ... ", GPU_tmp_b[i]);
          }
            CPUsum_of_GPU_tmp_b += GPU_tmp_b[i];
        }
        if (verbose > 1) printf("\n");
        printf("CPUsum_of_GPU_tmp_b = %f\n", CPUsum_of_GPU_tmp_b);

        free(GPU_tmp_b);
    }

  }

  // ------------ Final reduction: compute GPU_b with the reduction kernel ------------
  TIMER_START;

  {
      dim3 block_size(BLK_SIZE, 1, 1);
      dim3 grid_size(1, 1, 1);
      printf("%d: block_size = %d, grid_size = %d, kernel_len = %d\n", __LINE__, block_size.x, grid_size.x, kernel_len);
      reduction_kernel<<<grid_size, block_size>>>(kernel_len, dev_tmp_b, dev_b);
      bandwidth_numerator += (kernel_len+1)*sizeof(dtype);
  }


  checkCudaErrors( cudaDeviceSynchronize() );
  TIMER_STOP;
  *CompTime += TIMER_ELAPSED;
  *Bandwidth = bandwidth_numerator / ((*CompTime)*1e+9);

  // --------------- copy results from device to host ----------------

  dtype GPU_b;
  checkCudaErrors( cudaMemcpy(&GPU_b, dev_b, sizeof(dtype), cudaMemcpyDeviceToHost) );

  checkCudaErrors( cudaFree(dev_a) );
  checkCudaErrors( cudaFree(dev_b) );
  checkCudaErrors( cudaFree(dev_tmp_b) );
  return(GPU_b);
}


int main(int argc, char *argv[]) {

  printf("====================================== Problem computations ======================================\n");
// =========================================== Set-up the problem ============================================

  if (argc < 3) {
    printf("Usage: lab3_ex1 n v\n");
    return(1);
  }
  printf("argv[1] = %s\n", argv[1]);
  printf("argv[2] = %s\n", argv[2]);

  // ---------------- set-up the problem size -------------------

  int n = atoi(argv[1]), len = (1<<n), i;
  verbose = atoi(argv[2]);

  printf("n = %d --> len = 2^(n) = %d\n", n, len);
  printf("verbose = %d\n", verbose);
  printf("dtype = %s\n", XSTR(dtype));

//   // BUG check: code to generalize
//   if ((len%BLK_SIZE) != 0) {
//     printf("Now the code only support len (mod BLK_SIZE) == 0 (BLK_SIZE = %d)\n", BLK_SIZE);
//     exit(42);
//   }

  // ------------------ set-up the timers ---------------------

  TIMER_DEF;
  const char* lables[NPROBS] = {"CPU check", "kernel warp"};
  float errors[NPROBS], Times[NPROBS], Bandwidths[NPROBS], error;
  for (i=0; i<NPROBS; i++) {
    errors[i] = 1<<30;
    Bandwidths[i] = 0;
    Times[i] = 0;
  }


  // ------------------- set-up the problem -------------------

  dtype *a, GPU_b = 0.0f;
  double CPU_b = 0.0f;
  a = (dtype*)malloc(sizeof(dtype)*len);
  time_t t;
  srand((unsigned) time(&t));


#ifdef ONES_INIT
  for (i=0; i<len; i++)
    a[i] = 0.125f;
#endif

#ifdef DEBUG
  if (verbose > 1) {
    printf("a: ");
    for (i=0; i<len; i++)
        if (i==0 || a[i] != a[i-1])
          printf("%f ... ", a[i]);
    printf("\n");
  }
#endif
  // ======================================== Running the computations =========================================

  /* [ ... ]
   */

  // ========================== CPU computation =========================

  TIMER_START;
  for (i=0; i<len; i++) {
    CPU_b += a[i];
//     if (CPU_b == 16777216.0f) {
//       printf("ERROR: float overflow\n");
//       exit(42);
//     }
  }
  TIMER_STOP;

  errors[0] = 0.0f;
  Bandwidths[0] = 0.0f;
  Times[0] = TIMER_ELAPSED;

  printf("CPU_b = %f\n", CPU_b);

  printf("=========================== GPU Kernel Reduction Warp ===========================\n");
 
  GPU_b = execute_reduction_kernel(len, a, kernel_reduction_warp, &Bandwidths[1], &Times[1]);

  // ------------- Compare GPU and CPU solution --------------

  error = (float)fabs(CPU_b - GPU_b);
  errors[1] = error;

  printf("GPU_b = %f\n", GPU_b);

 

  // =========================== Print the performance results ===========================
  printf("Solution\n %9s\t%9s\t%9s\t%9s\n", "type", "error", "time (s)", "bandwidth (GB/s)");
  for (int i=0; i<NPROBS; i++) {
    printf("%12s:\t%9.7f\t%9.7f\t%9.7f\n", lables[i], errors[i], Times[i], Bandwidths[i]);
  }
  printf("\n");


  return(0);
}
