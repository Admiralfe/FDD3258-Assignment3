#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <mpi.h>

#define SEED     921
#define NUM_ITER 1000000000

int main(int argc, char* argv[])
{
    int provided;

    MPI_Init_thread(&argc, &argv, MPI_THREAD_SINGLE, &provided);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    printf("Rank %d says hello!\n", rank);

    int count = 0;
    double x, y, z, pi;
    srand(SEED * rank); // Important: Multiply SEED by "rank" when you introduce MPI!
    
    int iters_per_rank = NUM_ITER / size;
    int rank_0_iters = NUM_ITER - iters_per_rank * (size - 1);

    if (rank == 0) {
        for (int iter = 0; iter < rank_0_iters; iter++)
        {
            // Generate random (X,Y) points
            x = (double)random() / (double)RAND_MAX;
            y = (double)random() / (double)RAND_MAX;
            z = sqrt((x*x) + (y*y));
            
            // Check if point is in unit circle
            if (z <= 1.0)
            {
                count++;
            }
        }
    }
    
    else {
        // Calculate PI following a Monte Carlo method
        for (int iter = 0; iter < iters_per_rank; iter++)
        {
            // Generate random (X,Y) points
            x = (double)random() / (double)RAND_MAX;
            y = (double)random() / (double)RAND_MAX;
            z = sqrt((x*x) + (y*y));
            
            // Check if point is in unit circle
            if (z <= 1.0)
            {
                count++;
            }
        }
    }

    //Rank 0 collects results
    if (rank == 0) {
        for (int i = 1; i < size; ++i) {
            int tmp;
            MPI_Recv(&tmp, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            count += tmp;
        }
        // Estimate Pi and display the result
        pi = ((double)count / (double)NUM_ITER) * 4.0;
        printf("The result is %f\n", pi);
    }
    else {
        MPI_Send(&count, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
    }
    
    MPI_Finalize();
    return 0;
}