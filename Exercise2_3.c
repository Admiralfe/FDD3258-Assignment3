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

    int counts[size - 1];
    MPI_Request requests[size - 1];
    double start_time = MPI_Wtime();
    if (rank == 0) {
        //Post the non-blocking receive requests, so that receiving and computing can be done in parallell for rank 0...
        for (int sender = 1; sender < size; ++sender)
            MPI_Irecv(&counts[sender - 1], 1, MPI_INT, sender, 1, MPI_COMM_WORLD, &requests[sender - 1]);


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

        MPI_Waitall(size - 1, requests, MPI_STATUSES_IGNORE);
        //Rank 0 should have received all counts now
        for (int i = 0; i < size - 1; ++i)
            count += counts[i];
        
        pi = ((double)count / (double)NUM_ITER) * 4.0;
        printf("The result is %f\n", pi);

        double end_time = MPI_Wtime();
        printf("Time elapsed: %f s\n", end_time - start_time);
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
        MPI_Send(&count, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
    }
    

    MPI_Finalize(); 
    return 0;
}

