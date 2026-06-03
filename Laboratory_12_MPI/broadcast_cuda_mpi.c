#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <cuda_runtime.h>

// =========================
// Binomial tree broadcast
// =========================

void tree_bcast_gpu(int *d_buf, int count, int root, MPI_Comm comm) {
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    /* |========================================| */
    /* |           Put here your code           | */
    /* |========================================| */

}
// =========================
// Check correctness
// =========================
int check_buffer(int *buf, int count) {
    for (int i = 0; i < count; i++) {
        if (buf[i] != i) return 0;
    }
    return 1;
}
// =========================
// Main
// =========================
int main(int argc, char *argv[]) {
    int rank, size;
    int count = 1<<20;  // number of ints in the buffer
    int K = 50;       // repetitions for timing
    int root = 0;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); // global rank
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Local rank for GPU mapping
    int local_rank;
    MPI_Comm local_comm;

    MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0,
                        MPI_INFO_NULL, &local_comm);

    MPI_Comm_rank(local_comm, &local_rank);

    // GPU selection
    cudaSetDevice(local_rank);

    printf("Rank %d using GPU %d\n", rank, local_rank);

    if (argc == 3) {
        count = atoi(argv[1]);
        K     = atoi(argv[2]);
    } else if (rank == 0) {
        printf("Usage: mpirun -np P ./bcast_test count K\n");
        printf("Using defaults: count = %d, K = %d\n", count, K);
    }

    int *h_buf = (int *)malloc(count * sizeof(int));
    if (!h_buf) {
        fprintf(stderr, "Rank %d: cannot allocate buffer of size %d\n", rank, count);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    int *d_buf=NULL;

    cudaMalloc((void**)&d_buf, count * sizeof(int));

    /* =========================
       Tree broadcast (Isend/Irecv), staging as the baseline
       ========================= */
    // Initialize buffer on root and GPU
    if (rank == root) {
        for (int i = 0; i < count; i++) h_buf[i] = i;
    } else {
        for (int i = 0; i < count; i++) h_buf[i] = 0;
    }

    cudaMemcpy(d_buf, h_buf, count*sizeof(int), cudaMemcpyHostToDevice);
    cudaDeviceSynchronize();

    // Timing loop
    MPI_Barrier(MPI_COMM_WORLD);
    double t0 = MPI_Wtime();

    for (int k = 0; k < K; k++) {
        
        // GPU → CPU staging out 
        /* |========================================| */
        /* |           Put here your code           | */
        /* |========================================| */

        // Implement tree-based broadcast function using MPI_Isend, MPI_Irecv
        tree_bcast_gpu(/*Input the parameters*/); 

        // CPU → GPU staging in
        /* |========================================| */
        /* |           Put here your code           | */
        /* |========================================| */

    }

    MPI_Barrier(MPI_COMM_WORLD);
    double t1 = MPI_Wtime();
    double time_staging = (t1 - t0) / K;

    cudaMemcpy(h_buf, d_buf, count*sizeof(int), cudaMemcpyDeviceToHost);
    int ok_stage = check_buffer(h_buf, count);

    
    
    /* =========================
       EX3.a: Tree broadcast (Isend/Irecv), CUDA-aware
       ========================= */
    // Initialize buffer on root and GPU

    if (rank == root) {
        for (int i = 0; i < count; i++) h_buf[i] = i;
    } else {
        for (int i = 0; i < count; i++) h_buf[i] = 0;
    }
    
    cudaMemcpy(d_buf, h_buf, count * sizeof(int), cudaMemcpyHostToDevice);
    cudaDeviceSynchronize();

    // Timing loop
    MPI_Barrier(MPI_COMM_WORLD);
    t0 = MPI_Wtime();

    for (int k = 0; k < K; k++) {

        // Implement tree-based broadcast function using MPI_Isend, MPI_Irecv
        tree_bcast_gpu(/*Input the parameters*/); 

    }

    MPI_Barrier(MPI_COMM_WORLD);
    t1 = MPI_Wtime();
    double time_tree = (t1 - t0) / K;
   
    // Check correctness (once is enough)
    cudaMemcpy(h_buf, d_buf, count * sizeof(int), cudaMemcpyDeviceToHost);
    int ok_tree = check_buffer(h_buf, count);


    // =========================
    // EX3.b: MPI_Bcast (CUDA-aware)
    // =========================
    if (rank == root) {
        for (int i = 0; i < count; i++) h_buf[i] = i;
    } else {
        for (int i = 0; i < count; i++) h_buf[i] = 0;
    }

    cudaMemcpy(d_buf, h_buf, count*sizeof(int), cudaMemcpyHostToDevice);
    cudaDeviceSynchronize();

    MPI_Barrier(MPI_COMM_WORLD);
    t0 = MPI_Wtime();

    for (int k = 0; k < K; k++) {

        // Implement CUDA-aware collective broadcast using MPI_Bcast

        /* |========================================| */
        /* |           Put here your code           | */
        /* |========================================| */
    }

    MPI_Barrier(MPI_COMM_WORLD);
    t1 = MPI_Wtime();
    double time_bcast = (t1 - t0) / K;

    cudaMemcpy(h_buf, d_buf, count*sizeof(int), cudaMemcpyDeviceToHost);
    int ok_bcast = check_buffer(h_buf, count);


    // =========================
    // Results
    // =========================
    if (rank == root) {
        printf("\n===== RESULTS =====\n");
        printf("Count = %d ints (~%.2f MB), K = %d\n",
               count, count*sizeof(int)/1e6, K);

        printf("\nTree (P2P):        %e s  [%s]\n",
               time_tree, ok_tree ? "OK" : "WRONG");

        printf("MPI_Bcast (GPU):  %e s  [%s]\n",
               time_bcast, ok_bcast ? "OK" : "WRONG");

        printf("Staging (CPU):    %e s  [%s]\n",
               time_staging, ok_stage ? "OK" : "WRONG");

        printf("\nSpeedups:\n");
        printf("  Staging / Tree:   %.2fx\n", time_staging / time_tree);
        printf("  Staging / Bcast: %.2fx\n", time_staging / time_bcast);
    }

    free(h_buf);
    cudaFree(d_buf);
    MPI_Finalize();
    return 0;
}

