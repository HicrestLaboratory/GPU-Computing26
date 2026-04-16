# Laboratory 8

## Exercise 5

Given the access pattern represented in the figure implement a kernel that uses it to perform vector addition.

1. Write three different CUDA kernels, implement three memory access patterns;
2. Write an vector addition function on the CPU side, and measure its CPU_time. The purpose is to validate the correctness of GPU results;
3. Set 10 iterations for each kernel execution, measure their GPU_time and effective bandwidth;
4. Report the averages and standard deviation, compare the performance;


## Exercise 6
Find the best block and access pattern to perform an array element-wise addition kernel. Based on Exercise 5 (layout1 and layout 3), for each layout:
1. Fix the grid size to 16, vary the array size, block size:
    array size: 2^(12, 14, 16, 18, 20, 22)
    block size: 32, 64, 128, 256, 512

2. Report the average time solution (ms), std and bandwidth (GB/s)

3. Compare the performance, find the best block size and access pattern 








