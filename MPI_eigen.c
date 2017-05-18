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

/* Functions to normalize a vector */

double compute_norm(int vect_length,double x[vect_length]){
       
    int i;
    double norm, norm_sqr;
    norm = 0;
    norm_sqr = 0;    

    for (i = 0; i < vect_length; i++){
        norm_sqr += x[i]*x[i] ;
    }
    
    norm = sqrt(norm_sqr); 
    return norm;
}

void norm_vect(int vect_length,double x[vect_length]){
    
    int i;
    double  norm ;
    norm = compute_norm(vect_length,x) ;
   
    for (i = 0; i < vect_length; i++){
        x[i] = x[i] / norm ;
    }  
}


/* Function to multiply the sub row-matrix by the vector x and store it in x*/

void mat_mult(int heigth, int length, double A_rows[length*heigth], double x[length], double y[heigth]){

    int  i, j;
    
    for(i = 0; i < heigth; i++){
        for(j = 0; j < length ; j++){
           y[i] += A_rows[i*length + j] * x[j] ;
        }
    }    
}


int main(int argc, char *argv[]){

  int rank, nprocs, mat_size, blockrows_size, i, j, k, M ;
  double lambda, lambda0, e; 
 
  mat_size = atoi(argv[1]);
   
  //double x[2] = {1.0,1.0}; test of the norm_vect function


  double *A_blockrows;
  double *A;
  double *x;
  double *y;
  double *y_blocks;

  

  x = malloc(mat_size*sizeof(double));
  
  M = 2;
  e = 1e-6; //e is the tolerance
  lambda = 0;
  lambda0 = lambda + 2*e;
  k = 0;
  
  

  MPI_Request req_A_send;
  MPI_Request req_A_recv;

  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  blockrows_size = mat_size / nprocs ;


  y = malloc(mat_size*sizeof(double));
  y_blocks = malloc(blockrows_size*sizeof(double));

  /* Define new type to send datas */
  MPI_Datatype blocktype, blockselect;

  MPI_Type_contiguous(mat_size*blockrows_size, MPI_DOUBLE, &blocktype);
  MPI_Type_commit(&blocktype);

  MPI_Type_vector(blockrows_size, mat_size, mat_size, MPI_DOUBLE, &blockselect);
  MPI_Type_commit(&blockselect);

  if(rank == 0) {
    /*We initialize the matrix A and send rows of A to the other processes */
    A = (double *)malloc(sizeof(double) * mat_size * mat_size);
    for(i = 0; i < mat_size; i++){
        x[i]=1.0;
    }
    
    Prvalues(mat_size,1,x);
    norm_vect(mat_size,x);
    Prvalues(mat_size,1,x);

    
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
     
     /* double norm = compute_norm(2,x);
      printf("norm %.3f\t \n",norm);         // Testing the normalization of the vectors
      norm_vect(2,x);
      Prvalues(1,2,x);*/       
  }

  

  /* After receiving rows of A, we start the power method and use Allgather to send y to all processes */
  A_blockrows = (double *)malloc(mat_size * blockrows_size * sizeof(double));
  MPI_Irecv(A_blockrows, 1, blocktype, 0, 3, MPI_COMM_WORLD, &req_A_recv);

  MPI_Wait(&req_A_recv, MPI_STATUS_IGNORE);
  

  MPI_Bcast(&x[0],mat_size,MPI_DOUBLE,0,MPI_COMM_WORLD);
  //if(rank == 2)mat_mult( blockrows_size,mat_size,A_blockrows,x,y);

  
  while (abs(lambda-lambda0) >= e && k <= M ){

    k += 1;
  }
  
  
  

  /* Printing Part*/
  sleep(rank);
  printf("rank %d \n", rank);
  Prvalues(mat_size,1,x);
  //if (rank == 2) Prvalues(1, blockrows_size, y);
  Prvalues(mat_size, blockrows_size, A_blockrows);

  /*When the algorithm terminates we print the dominant eigein value and the associated eigen vector in root processor */ 

  MPI_Type_free(&blockselect);
  MPI_Type_free(&blocktype);
  
  MPI_Finalize();


  return 0;

}
