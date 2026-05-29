#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <cuda_runtime.h>

void linear_bcast_blocking_gpu(int *h_buf, int count, int root, MPI_Comm comm) {
    /* |========================================| */
    /* |           Put here your code           | */
    /* |========================================| */
}

void tree_bcast_nonblocking_gpu(int *h_buf, int count, int root, MPI_Comm comm) {
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    /* |========================================| */
    /* |           Put here your code           | */
    /* |========================================| */
}

int check_buffer(int *buf, int count) {
    for (int i = 0; i < count; i++) {
        if (buf[i] != i) return 0;
    }
    return 1;
}

int main(int argc, char *argv[]) {
    int rank, size;
    int count = 1024;  // number of ints in the buffer
    int K = 100;       // repetitions for timing
    int root = 0;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); // global rank
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Map 4 ranks to 4 GPUs
    int local_rank;
    MPI_Comm local_comm;

    /* |========================================| */
    /* |           Put here your code           | */
    /* |========================================| */

    printf("Rank %d using GPU %d\n", rank, local_rank);

    // Initialize buffer
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
       Ex1.a: Linear broadcast (blocking)
       ========================= */
    // Initialize buffer on root and GPU
    if (rank == root) {
        for (int i = 0; i < count; i++) h_buf[i] = i;
    }
    cudaMemcpy(d_buf, h_buf, count * sizeof(int), cudaMemcpyHostToDevice);

    // Ex1.a: Timing loop
    MPI_Barrier(MPI_COMM_WORLD);
    double t0 = MPI_Wtime();

    for (int k = 0; k < K; k++) {
        // GPU → CPU staging out 
        /* |========================================| */
        /* |           Put here your code           | */
        /* |========================================| */
        // linear broadcast (blocking) function
        linear_bcast_blocking_gpu(h_buf, count, root, MPI_COMM_WORLD);
        
        // CPU → GPU staging in
        /* |========================================| */
        /* |           Put here your code           | */
        /* |========================================| */
    }

    MPI_Barrier(MPI_COMM_WORLD);
    double t1 = MPI_Wtime();
    double time_linear = (t1 - t0) / K;

    // Check correctness (once is enough)
    cudaMemcpy(h_buf, d_buf, count * sizeof(int), cudaMemcpyDeviceToHost);
    int ok_linear = check_buffer(h_buf, count);

    if (!ok_linear) {
        printf("Rank %d: linear broadcast produced WRONG data!\n", rank);
    }

    if (rank == root) {
        printf("Linear broadcast: count = %d, K = %d\n", count, K);
        printf("  Avg time per broadcast: %e seconds\n", time_linear);
    }

    /* =========================
       Ex1.b:Tree broadcast (Isend/Irecv)
       ========================= */
    // Reinitialize buffer on root and GPU

    if (rank == root) {
        for (int i = 0; i < count; i++) h_buf[i] = i;
    } else {
        for (int i = 0; i < count; i++) h_buf[i] = 0;
    }
    cudaMemset(d_buf, 0, count*sizeof(int));
    cudaMemcpy(d_buf, h_buf, count * sizeof(int), cudaMemcpyHostToDevice);

    // Ex1.b:Timing loop
    MPI_Barrier(MPI_COMM_WORLD);
    double t2 = MPI_Wtime();

    for (int k = 0; k < K; k++) {

        // GPU → CPU staging out 
        /* |========================================| */
        /* |           Put here your code           | */
        /* |========================================| */

        tree_bcast_nonblocking_gpu(h_buf, count, root, MPI_COMM_WORLD);

        // CPU → GPU staging in
        /* |========================================| */
        /* |           Put here your code           | */
        /* |========================================| */
    }

    MPI_Barrier(MPI_COMM_WORLD);
    double t3 = MPI_Wtime();
    double time_tree = (t3 - t2) / K;
   
    // Check correctness (once is enough)
    cudaMemcpy(h_buf, d_buf, count * sizeof(int), cudaMemcpyDeviceToHost);
    int ok_tree = check_buffer(h_buf, count);
    if (!ok_tree) {
        printf("Rank %d: tree broadcast produced WRONG data!\n", rank);
    }

    if (rank == root) {
        printf("Tree (Isend/Irecv) broadcast: count = %d, K = %d\n", count, K);
        printf("  Avg time per broadcast: %e seconds\n", time_tree);
        printf("Speedup (linear / tree): %.2fx\n", time_linear / time_tree);
    }

    free(h_buf);
    cudaFree(d_buf);
    MPI_Finalize();
    return 0;
}

