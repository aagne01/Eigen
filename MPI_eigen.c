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

 
  }
  
  MPI_Finalize();


  return 0;

}
