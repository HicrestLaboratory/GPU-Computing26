#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>

int get_random() {
    return rand() % 100;
}

int main(int argc, char** argv) {
    // Initialize MPI
    MPI_Init(&argc, &argv);
    int my_rank, comm_size; 
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
    int number = get_random();
    // Send Recv communication between two processes
    /* |========================================| */
    /* |           Put here your code           | */
    /* |========================================| */
    
    MPI_Finalize();
    return 0;
}
