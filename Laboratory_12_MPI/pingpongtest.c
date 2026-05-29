#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <cuda_runtime.h>

int main(int argc, char *argv[]) {
    int rank, size;
    int N = 10;   // iterations
    int h_msg;
    int *d_msg = NULL;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);  // Get the rank of the process
    MPI_Comm_size(MPI_COMM_WORLD, &size);// Get total number of processes

    if (rank == 0) {
        printf("Running CPU <-> GPU ping-pong with 1 MPI process\n");
    

        // Select GPU (important!)
        int gpu_id = 0;  // you can change this
        cudaSetDevice(gpu_id);

        // Allocate GPU memory
        cudaMalloc((void**)&d_msg, sizeof(int));

        for (int i = 0; i < N; i++) {
            h_msg = i;

            // CPU → GPU
            cudaMemcpy(d_msg, &h_msg, sizeof(int), cudaMemcpyHostToDevice);

            // GPU → CPU
            cudaMemcpy(&h_msg, d_msg, sizeof(int), cudaMemcpyDeviceToHost);

            printf("Iteration %d, value = %d\n", i, h_msg);
        }
        cudaFree(d_msg);
    }
    
    MPI_Finalize();
    return 0;
}