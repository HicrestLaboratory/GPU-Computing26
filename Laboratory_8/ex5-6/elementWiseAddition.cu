#include <stdio.h>
#include <time.h>
#include <sys/time.h>
#include <math.h>
#include <cuda_runtime.h>

#define NITER 10
#define WARMUP 4

#define NPROBS 3

#define TIMER_DEF     struct timeval temp_1, temp_2

#define TIMER_START   gettimeofday(&temp_1, (struct timezone*)0)

#define TIMER_STOP    gettimeofday(&temp_2, (struct timezone*)0)

#define TIMER_ELAPSED ((temp_2.tv_sec-temp_1.tv_sec)*1000.0 +(temp_2.tv_usec-temp_1.tv_usec)/1000.0) // return ms

#define DBG_CHECK { printf("DBG_CHECK: file %s at line %d\n", __FILE__, __LINE__ ); }
// #define DEBUG
// #define BLK_DISPACH

#define RUN_SOLUTIONS

#define BLK_SIZE 32
#define GRD_SIZE 16

#define STR(s) #s
#define XSTR(s) STR(s)
#define dtype float

// Print device properties
void printDevProp(cudaDeviceProp devProp)
{
    printf("Major revision number:         %d\n",  devProp.major);
    printf("Minor revision number:         %d\n",  devProp.minor);
    printf("Name:                          %s\n",  devProp.name);
    printf("  Memory Clock rate:           %.0f Mhz\n", devProp.memoryClockRate * 1e-3f);

    printf("  Memory Bus Width:            %d bit\n",devProp.memoryBusWidth);

    printf("  Peak Memory Bandwidth:       %7.3f GB/s\n",2.0*devProp.memoryClockRate*(devProp.memoryBusWidth/8)/1.0e6);

    printf("  Multiprocessors:             %3d\n",devProp.multiProcessorCount);
    printf("  Maximum number of threads per multiprocessor:  %d\n",devProp.maxThreadsPerMultiProcessor);
    printf("  Maximum number of threads per block:           %d\n",devProp.maxThreadsPerBlock);
    printf("  Max dimension size of a thread block (x,y,z): (%d, %d, %d)\n",
           devProp.maxThreadsDim[0], devProp.maxThreadsDim[1],devProp.maxThreadsDim[2]);
    printf("  Max dimension size of a grid size    (x,y,z): (%d, %d, %d)\n",
           devProp.maxGridSize[0], devProp.maxGridSize[1],devProp.maxGridSize[2]);
    printf("  Total amount of shared memory per block:       %zu bytes\n", devProp.sharedMemPerBlock);
    return;
}

float arithmetic_mean(float *v, int len) {

    /* |========================================| */
    /* |           Put here your code           | */
    /* |========================================| */

}

float compute_std(float *arr, int n, float mean) {
    /* |========================================| */
    /* |           Put here your code           | */
    /* |========================================| */

}


__device__ uint get_smid(void) {

     uint ret;

     asm("mov.u32 %0, %smid;" : "=r"(ret) );

     return ret;

}

__global__
void example_kernel(int n, dtype *a, dtype* b, dtype* c)
{
  if (threadIdx.x==0)
      printf("block %d runs on sm %d\n", blockIdx.x, get_smid());

  // [ ... ]
}



__global__
void Problem_solution_layout1(int n, dtype *a, dtype* b, dtype* c)
{
    /* |========================================| */
    /* |           Put here your code           | */
    /* |========================================| */

}

__global__
void Problem_solution_layout2(int n, dtype *a, dtype* b, dtype* c, int h)
{
    /* |========================================| */
    /* |           Put here your code           | */
    /* |========================================| */

}

__global__
void Problem_solution_layout3(int n, dtype *a, dtype* b, dtype* c)
{
    /* |========================================| */
    /* |           Put here your code           | */
    /* |========================================| */

}


int main(int argc, char *argv[]) {

  printf("======================================= Device properties ========================================\n");

  int deviceCount = 0;
  cudaError_t error_id = cudaGetDeviceCount(&deviceCount);

  int dev, SM_num;
  for (dev = 0; dev < deviceCount; ++dev) {
    cudaSetDevice(dev);
    cudaDeviceProp deviceProp;
    cudaGetDeviceProperties(&deviceProp, dev);
    printDevProp(deviceProp);
    SM_num = deviceProp.multiProcessorCount;

  }


  printf("====================================== Problem computations ======================================\n");
// =========================================== Set-up the problem ============================================

  if (argc < 2) {
    printf("Usage: lab8_ex5 n\n");
    return(1);
  }
  printf("argv[1] = %s\n", argv[1]);

  // ---------------- set-up the problem size -------------------

  int n = atoi(argv[1]), len = (1<<n), i;

  printf("n = %d --> len = 2^(n) = %d\n", n, len);
  printf("dtype = %s\n", XSTR(dtype));


  // ------------------ set-up the variables ---------------------

  TIMER_DEF;
  float error, cputime, gputime_iter; 
  float gputime_layout1[NITER], gputime_layout2[NITER], gputime_layout3[NITER];
  float errors[3], gputimes[3], std[3], bandwidth[3];

  // ------------------- set-up the problem -------------------

  dtype *a, *b, *CPU_c, *GPU_c;
  a = (dtype*)malloc(sizeof(dtype)*len);
  b = (dtype*)malloc(sizeof(dtype)*len);
  CPU_c = (dtype*)malloc(sizeof(dtype)*len);
  GPU_c = (dtype*)malloc(sizeof(dtype)*len);
  time_t t;
  srand((unsigned) time(&t));

  int typ = (strcmp( XSTR(dtype) ,"int")==0);
  if (typ) {
      // here we generate random ints
      int rand_range = (1<<11);
      printf("rand_range= %d\n", rand_range);
      for (i=0; i<len; i++) {
          a[i] = rand()/(rand_range);
          b[i] = rand()/(rand_range);
          GPU_c[i] = (dtype)0;
      }
  } else {
      // here we generate random floats
      for (i=0; i<len; i++) {
        a[i] = (dtype)rand()/((dtype)RAND_MAX);
        b[i] = (dtype)rand()/((dtype)RAND_MAX);
        GPU_c[i] = (dtype)0;
      }
  }


// ======================================== EX5: Implement three access patterns =========================================

  // ---------------- computation on CPU -------------------
  // 1) Write an vector addition function on the CPU side, and measure its CPU_time
    
    /* |========================================| */
    /* |           Put here your code           | */
    /* |========================================| */





  // ---------------- allocating GPU vectors -------------------
  dtype *dev_a, *dev_b, *dev_c;

  cudaMalloc(&dev_a, len*sizeof(dtype));
  cudaMalloc(&dev_b, len*sizeof(dtype));
  cudaMalloc(&dev_c, len*sizeof(dtype));

    // ------------ copy date from host to device --------------

  cudaMemcpy(dev_a, a, len*sizeof(dtype), cudaMemcpyHostToDevice);
  cudaMemcpy(dev_b, b, len*sizeof(dtype), cudaMemcpyHostToDevice);

  dim3 block_size(BLK_SIZE, 1, 1);
  dim3 grid_size(GRD_SIZE, 1, 1);
  printf("block_size = %d, grid_size = %d, elements per thread = %f\n", block_size.x, grid_size.x, (float)len/(block_size.x*grid_size.x));
  
  
  // Layout1 implementation: Coalesced access, including the following three steps
  
  // Step1: ------------ computation solution with Layout 1 -----------

    /* |========================================| */
    /* |           Put here your code           | */
    /* |========================================| */

  
  // Step2: ----------- copy results from device to host ------------

    /* |========================================| */
    /* |           Put here your code           | */
    /* |========================================| */

  // Step3: ------------- Compare GPU and CPU solution, Report the time/std/bandwidth --------------

    /* |========================================| */
    /* |           Put here your code           | */
    /* |========================================| */


  // ------------------- reset gpu buffers ---------------------

  cudaMemset(dev_a, 0, len*sizeof(dtype));
  cudaMemset(dev_b, 0, len*sizeof(dtype));
  cudaMemcpy(dev_a, a, len*sizeof(dtype), cudaMemcpyHostToDevice);
  cudaMemcpy(dev_b, b, len*sizeof(dtype), cudaMemcpyHostToDevice);


  // Layout2 implementation: Coalesced access + offset,

  // Step1: ------------ computation solution with Layout 2 -----------

  /* |========================================| */
  /* |           Put here your code           | */
  /* |========================================| */


// Step2: ----------- copy results from device to host ------------

  /* |========================================| */
  /* |           Put here your code           | */
  /* |========================================| */

// Step3: ------------- Compare GPU and CPU solution, Report the time/std/bandwidth --------------

  /* |========================================| */
  /* |           Put here your code           | */
  /* |========================================| */


    // ------------------- reset gpu buffers ---------------------

  cudaMemset(dev_a, 0, len*sizeof(dtype));
  cudaMemset(dev_b, 0, len*sizeof(dtype));
  cudaMemcpy(dev_a, a, len*sizeof(dtype), cudaMemcpyHostToDevice);
  cudaMemcpy(dev_b, b, len*sizeof(dtype), cudaMemcpyHostToDevice);


  // Layout3 implementation: Stride access,

  // Step1: ------------ computation solution with Layout 3 -----------

  /* |========================================| */
  /* |           Put here your code           | */
  /* |========================================| */


// Step2: ----------- copy results from device to host ------------

  /* |========================================| */
  /* |           Put here your code           | */
  /* |========================================| */

// Step3: ------------- Compare GPU and CPU solution, Report the time/std/bandwidth --------------

  /* |========================================| */
  /* |           Put here your code           | */
  /* |========================================| */


// ============================================ Print the results ============================================
  printf("================================== Times and results of my code ==================================\n");
  printf("\nVector len = %d, CPU time(ms) = %5.3f\n\n", len, cputime);
  printf("\t\tError\tTime(ms)\tstd\tBandwidth(GB/s)\n");
  for (int i=0; i<3; i++)
    printf("Layout %d:\t%5.3f\t%5.3f\t%5.3f\t%5.3f\n", i+1, errors[i], gputimes[i],std[i],bandwidth[i]);

// ======================================== EX5: End =========================================

// ===================================== EX6: Find best BLK dimensions ====================================

  printf("===================================== Find best BLK dimensions ====================================\n");

  int blk_size; // the range of block size: 32, 64, 128, 256, 512
  float compareDimStd[5], compareDimTime[5],compareDimBw[5]; // to store the result of each blk_size

// Implement layout1 and layout3 separately

  /* |========================================| */
  /* |           Put here your code           | */
  /* |========================================| */



  printf("\nVector len = %d, CPU time(ms) = %5.5f, kernel = %s\n\n\n", len, cputime, kernel_string);

  printf("Block size\tTime(Std)\tBandwidth(GB/s)\n");
  for (int i=0; i<5; i++) {
    printf("%d:\t", 1 << (i+5));

    printf("%.3f(%4.2f)\t%.3f\t", compareDimTime[i], compareDimStd[i],compareDimBw[i]);
    printf("\n");
  }
  printf("\n");


  // ----------------- free GPU variable ---------------------

  cudaFree(dev_a);
  cudaFree(dev_b);
  cudaFree(dev_c);

  // ---------------------------------------------------------

  return(0);
}
