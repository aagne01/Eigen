#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <time.h>
#include <string.h>


int main(int argc, char *argv[]){

  int rank, nprocs ;
  

  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  if(rank == 0) {
  /*We initialize the matrix A and send rows of A to the other processes */
 
  }

  /* After receiving rows of A, we start the power method and use Allgather to send y to all processes */

  /*When the algorithm terminates we print the dominant eigein value and the associated eigen vector in root processor */ 
  
  MPI_Finalize();


  return 0;

}
