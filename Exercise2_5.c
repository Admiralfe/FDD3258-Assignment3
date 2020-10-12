#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <time.h>

#define SEED     921
#define NUM_ITER 1000000000

int main(int argc, char* argv[])
{
    int count = 0;
    double x, y, z, pi;

    int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_SINGLE, &provided);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    srand(rank * SEED); // Important: Multiply SEED by "rank" when you introduce MPI!

    int iters_per_rank = NUM_ITER / size;
    int rank_0_iters = NUM_ITER - iters_per_rank * (size - 1);

    int counts[size];
    double start_time = MPI_Wtime();
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

    int global_count;
    MPI_Reduce(&count, &global_count, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        pi = ((double)global_count / (double)NUM_ITER) * 4.0;
        printf("The result is %f\n", pi);

        double end_time = MPI_Wtime();
        printf("Time elapsed: %f s\n", end_time - start_time);
        
    }
        

    MPI_Finalize(); 
    return 0;
}

