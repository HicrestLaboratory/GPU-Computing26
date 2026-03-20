# Laboratory 3

## Exercise 1

Given a command-line parameter n, generate two random square matrices, perform a general matrix multiplication using OpenBLAS. 
1)Benchmark its runtime by performing two warm-up cycles and averaging over 10 rounds. 
2)Compute the average with both arithmetic and geometric mean.
3)Complete the provided Makefile and use the provided sbatch script.

You can find a template code inside "ex1/". Keep all the benchmark function into the library "ex1/src/my_time_lib.c" and "ex1/include/my_time_lib.h". 

## Exercise 2
Benchmark a program that use OpenBLAS to perform a general matrix multiplication between two random square matrices. Get the execution time and the FLOP/s.

Implement your own general matrix multiplication and compare your performance with the OpenBLAS implementation.

The template code is inside "ex2/"

## Exercise 3
Implement a simple GEMM algorithm that multiplies two matrices of size NxN by splitting the orginal matrix into blocks
1)The block_size is such that  NxN is a multiple of block_size X block_size
2)Measure the FLOPS of your sequential implementations
Note: Perform a correctness check. And for debugging, compile the code using –g option;

Run the code with valgrind to measure the cache hit/miss




