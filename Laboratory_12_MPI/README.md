# Laboratory 12

## Exercise 1: P2P Broadcast 

In this exercise, you will implement your own broadcast operation using point-to-point communication, and involve GPU buffer:

1. Set 4 MPI processes and map to 4 GPUs
2. A naive broadcast using blocking MPI_Send / MPI_Recv + staging
3. An optimized tree-based broadcast using non-blocking MPI_Isend / MPI_Irecv + MPI_Wait + staging.

Measure and print, for each implementation:
1. the total time for K broadcasts,
2. the average time per broadcast.
3. Verify correctness by having non-root processes check that they received the same data as the root.


## Exercise 2: Collective Broadcast

In this exercise, you will implement your own broadcast operation using collective communication, and involve GPU buffer:

1. Set 4 MPI processes and map to 4 GPUs
2. Use blocking MPI_Bcast + staging
3. Use non-blocking MPI_Ibcast+ MPI_Wait + staging
4. Compare their performance by measureing 









