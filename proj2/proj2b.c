#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "mpi.h"

//set the number of steps
#define NUMSTEPS 1000000

int main(int argc, char** argv) {
    int rank, size, i;
    double x, pi, sum = 0.0, total_sum = 0.0;
    //start timer
    struct timespec start, end;
    clock_gettime(CLOCK_MONOTONIC, &start);
  
    //Intialize the MPI
    MPI_Init(&argc, &argv);
    //get the rank
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    //get the total number of processes
    MPI_Comm_size(MPI_COMM_WORLD, &size);
  
    //create the step size
    double step = 1.0 / (double) NUMSTEPS; 
    //loop throught the local range
    for (i = rank; i < NUMSTEPS; i += size) { 
        x = (i + 0.5) * step;
        //calculate the local sum of the series
        sum += 4.0 / (1.0 + x * x); 
    }
    //check for root process
    if (rank == 0) {
        for (i = 1; i < size; i++) {
            //call for the total sum of other processes
            MPI_Recv(&total_sum, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            //add local sums to the final sum
            sum += total_sum;
        }
      //if not the root process
    } else {
        //send the local sum to the root process
        MPI_Send(&sum, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    }
    //check if the root process
    if (rank == 0) {
        //calculate the value of pi
        pi = step * sum;
        printf("PI is %.20f\n", pi);
        //end timer
        clock_gettime(CLOCK_MONOTONIC, &end);
        u_int64_t diff = 1000000000L * (end.tv_sec - start.tv_sec) + end.tv_nsec - start.tv_nsec;
        //display
        printf("elapsed time = %llu nanoseconds\n", (long long unsigned int) diff);
    }
    //close MPI
    MPI_Finalize();
    return 0;
}
