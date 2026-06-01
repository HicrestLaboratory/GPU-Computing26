#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <cuda_runtime.h>

void linear_bcast_blocking_gpu(int *h_buf, int count, int root, MPI_Comm comm) {
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    if (rank == root) {

        for (int dest = 0; dest < size; dest++) {
            if (dest == root) continue;
            MPI_Send(h_buf, count, MPI_INT, dest, 0, comm);
        }

    } else {
        MPI_Recv(h_buf, count, MPI_INT, root, 0, comm, MPI_STATUS_IGNORE);
    }
}

void tree_bcast_nonblocking_gpu(int *h_buf, int count, int root, MPI_Comm comm) {
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    int relative_rank = (rank - root + size) % size;

    MPI_Request req;
    int mask = 1;

    while (mask < size) {

        if (relative_rank < mask) { //senders
            int dst = relative_rank + mask;
            if (dst < size) {
                int dst_rank = (dst + root) % size;

                MPI_Isend(h_buf, count, MPI_INT, dst_rank, 0, comm, &req);
                MPI_Wait(&req, MPI_STATUS_IGNORE);
            }
        }
        else if (relative_rank < 2 * mask) { //receivers
            int src = relative_rank - mask;
            int src_rank = (src + root) % size;

            MPI_Irecv(h_buf, count, MPI_INT, src_rank, 0, comm, &req);
            MPI_Wait(&req, MPI_STATUS_IGNORE);
        }

        mask <<= 1;
    }
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
       Linear broadcast (blocking)
       ========================= */
    // Initialize buffer on root and GPU
    if (rank == root) {
        for (int i = 0; i < count; i++) h_buf[i] = i;
    }
    cudaMemcpy(d_buf, h_buf, count * sizeof(int), cudaMemcpyHostToDevice);

    // Timing loop
    MPI_Barrier(MPI_COMM_WORLD);
    double t0 = MPI_Wtime();

    for (int k = 0; k < K; k++) {

        if (rank == root) {
            cudaMemcpy(h_buf, d_buf, count * sizeof(int), cudaMemcpyDeviceToHost);   // GPU → CPU staging
        }

        linear_bcast_blocking_gpu(h_buf, count, root, MPI_COMM_WORLD);

        cudaMemcpy(d_buf, h_buf, count * sizeof(int), cudaMemcpyHostToDevice);      // CPU → GPU restore
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
       Tree broadcast (Isend/Irecv)
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

        cudaMemcpy(h_buf, d_buf, count * sizeof(int), cudaMemcpyDeviceToHost);   // stage out

        tree_bcast_nonblocking_gpu(h_buf, count, root, MPI_COMM_WORLD);

        cudaMemcpy(d_buf, h_buf, count * sizeof(int), cudaMemcpyHostToDevice);   // stage in
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

