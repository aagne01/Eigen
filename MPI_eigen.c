#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <time.h>
#include <string.h>



/* Function to Print matrix */
void Prvalues(int length, int heigth,  double matrix[length * heigth]){   
    int i, j;
    printf("\n");
    for (i = 0; i < heigth; i++){
        for (j = 0; j < length; j++){
            printf("%.1f\t", matrix[i*length + j]);
        }
        printf("\n");
    }
    printf("\n");
}



int main(int argc, char *argv[]){

  int rank, nprocs, mat_size, blockrows_size, i, j ;

  mat_size = atoi(argv[1]);


  double *A_blockrows;
  double *A = NULL;

  MPI_Request req_A_send;
  MPI_Request req_A_recv;

  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  blockrows_size = mat_size / nprocs ;

  /* Define new type to send datas */
  MPI_Datatype blocktype, blockselect;

  MPI_Type_contiguous(mat_size*blockrows_size, MPI_DOUBLE, &blocktype);
  MPI_Type_commit(&blocktype);

  MPI_Type_vector(blockrows_size, mat_size, mat_size, MPI_DOUBLE, &blockselect);
  MPI_Type_commit(&blockselect);

  if(rank == 0) {
    /*We initialize the matrix A and send rows of A to the other processes */
    A = (double *)malloc(sizeof(double) * mat_size * mat_size);

    
    /* Matrix generation*/
    for(i = 0; i < mat_size; i++){
       for(j = 0; j < mat_size; j++){
          A[i * mat_size + j] = i + 0.1*j;   
       }    
    }
    Prvalues(mat_size, mat_size, A);
 
    /* Distribution of Block_rows of A */
     for (i = 0; i < nprocs; i++){
            int i_b_m = i * blockrows_size * mat_size;
            MPI_Isend(&A[i_b_m], 1, blockselect, i, 3, MPI_COMM_WORLD, &req_A_send);
     }

  }


  /* After receiving rows of A, we start the power method and use Allgather to send y to all processes */
   A_blockrows = (double *)malloc(mat_size * blockrows_size * sizeof(double));
   MPI_Irecv(A_blockrows, 1, blocktype, 0, 3, MPI_COMM_WORLD, &req_A_recv);

   MPI_Wait(&req_A_recv, MPI_STATUS_IGNORE);


  /* Printing Part*/
  sleep(rank);
  printf("rank %d \n", rank);
  Prvalues(mat_size, blockrows_size, A_blockrows);

  /*When the algorithm terminates we print the dominant eigein value and the associated eigen vector in root processor */ 

  MPI_Type_free(&blockselect);
  MPI_Type_free(&blocktype);
  
  MPI_Finalize();


  return 0;

}
