#include <mpi.h>
#include <stdio.h>

int main(int argc, char** argv) {

    MPI_Init(&argc, &argv);

    int rank, size;
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    // Get the rank of the process
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    // Get total number of processes
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    // Get the name of the processor
    MPI_Get_processor_name(processor_name, &name_len); 

    printf("Hello from rank %d of %d running on %s\n",
           rank, size, processor_name);

    MPI_Finalize();
    return 0;
}

