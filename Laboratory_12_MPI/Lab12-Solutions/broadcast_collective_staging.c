#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <cuda_runtime.h>

void bcast_collective_blocking_gpu(int *h_buf, int *d_buf, int count, int root, MPI_Comm comm) {
    int rank;
    MPI_Comm_rank(comm, &rank);

    if (rank == root) {
        // GPU → host staging
        cudaMemcpy(h_buf, d_buf, count * sizeof(int), cudaMemcpyDeviceToHost);
    }

    // Collective communication
    MPI_Bcast(h_buf, count, MPI_INT, root, comm);

    // host → GPU staging
    cudaMemcpy(d_buf, h_buf, count * sizeof(int), cudaMemcpyHostToDevice);
}

void bcast_collective_nonblocking_gpu(int *h_buf, int *d_buf, int count,
                                      int root, MPI_Comm comm, cudaStream_t stream) {
    int rank;
    MPI_Comm_rank(comm, &rank);

    MPI_Request req;

    if (rank == root) {
        cudaMemcpyAsync(h_buf, d_buf, count * sizeof(int), cudaMemcpyDeviceToHost, stream);

        cudaStreamSynchronize(stream);
    }

    MPI_Ibcast(h_buf, count, MPI_INT, root, comm, &req);

    MPI_Wait(&req, MPI_STATUS_IGNORE);

    cudaMemcpyAsync(d_buf, h_buf, count * sizeof(int), cudaMemcpyHostToDevice, stream);

    cudaStreamSynchronize(stream);
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
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
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
    int *d_buf = NULL;

    cudaMalloc((void**)&d_buf, count * sizeof(int));


    /* =========================
       Collective broadcast (blocking)
       ========================= */
    // Initialize buffer on root and GPU
    if (rank == root) {
        for (int i = 0; i < count; i++) h_buf[i] = i;
    } else {
        for (int i = 0; i < count; i++) h_buf[i] = 0;
    }
    cudaMemcpy(d_buf, h_buf, count * sizeof(int), cudaMemcpyHostToDevice);

    // Timing loop
    cudaStream_t stream;
    cudaStreamCreate(&stream);

    MPI_Barrier(MPI_COMM_WORLD);
    double t0 = MPI_Wtime();

    for (int k = 0; k < K; k++) {
        bcast_collective_blocking_gpu(h_buf, d_buf, count, root, MPI_COMM_WORLD);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    double t1 = MPI_Wtime();
    double time_blocking = (t1 - t0) / K;

    // Check correctness (once is enough)
    cudaMemcpy(h_buf, d_buf, count * sizeof(int), cudaMemcpyDeviceToHost);
    
    int ok_blocking = check_buffer(h_buf, count);
    if (!ok_blocking) {
        printf("Rank %d: collective broadcast blocking produced WRONG data!\n", rank);
    }

    if (rank == root) {
        printf("Collective broadcast blocking: count = %d, K = %d\n", count, K);
        printf("  Avg time per broadcast: %e seconds\n", time_blocking);
    }

    /* =========================
       Collective broadcast (non-blocking)
       ========================= */
    // Reinitialize buffer on root and GPU

    if (rank == root) {
        for (int i = 0; i < count; i++) h_buf[i] = i;
    } else {
        for (int i = 0; i < count; i++) h_buf[i] = 0;
    }
    cudaMemset(d_buf, 0, count*sizeof(int));
    cudaMemcpy(d_buf, h_buf, count * sizeof(int), cudaMemcpyHostToDevice);

    // Timing loop
    MPI_Barrier(MPI_COMM_WORLD);
    double t2 = MPI_Wtime();

    for (int k = 0; k < K; k++) {
        bcast_collective_nonblocking_gpu(h_buf, d_buf, count, root,
                                        MPI_COMM_WORLD, stream);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    double t3 = MPI_Wtime();
    double time_nonblocking = (t3 - t2) / K;

    // Check correctness (once is enough)
    cudaMemcpy(h_buf, d_buf, count * sizeof(int), cudaMemcpyDeviceToHost);
    
    int ok_nonblocking = check_buffer(h_buf, count);
    if (!ok_nonblocking) {
        printf("Rank %d: collective broadcast non-blocking produced WRONG data!\n", rank);
    }

    if (rank == root) {
        printf("Collective broadcast non-blocking: count = %d, K = %d\n", count, K);
        printf("  Avg time per broadcast: %e seconds\n", time_nonblocking);
        printf("Speedup (blocking / nonblocking): %.2fx\n", time_blocking / time_nonblocking);
    }

    free(h_buf);
    cudaFree(d_buf);
    cudaStreamDestroy(stream);
    MPI_Finalize();
    return 0;
}

