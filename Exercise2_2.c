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

    //Binary tree reduction
    for (int step = 2; step <= size; step *= 2) {
        if ((rank % step) == 0) {
            int other_count;
            MPI_Recv(&other_count, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            count += other_count;
        }
        
        if ((rank % step) == step / 2) {
            MPI_Send(&count, 1, MPI_INT, rank - (step / 2), 0, MPI_COMM_WORLD);
        }
    }
    
    //Rank 0 should now have final counts
    if (rank == 0) {
        pi = ((double) count / (double) NUM_ITER) * 4.0;
        printf("The result is %f\n", pi);
        double end_time = MPI_Wtime();
        printf("Execution time: %f s\n", end_time - start_time);
    }

    MPI_Finalize();
    return 0;
}