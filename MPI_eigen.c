#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <time.h>
#include <string.h>


void Prvalues(int length, int heigth,  double matrix[length * heigth]);
double compute_norm(int vect_length,double x[vect_length]);
void norm_vect(int vect_length,double x[vect_length]);
void mat_mult(int length, double A_rows[], double x[], double y[]);
double power_method(double A[],int n_iterations,int mat_size);
void generate_matrix(double matrix[],int mat_size);


int main(int argc, char *argv[]){

  int rank, nprocs, mat_size, blockrows_size, i, j, M ;
  //double lambda, e;
 
  mat_size = atoi(argv[1]);
   

/*  double *A_blockrows;
  double *A;
 */ 
  
  M = atoi(argv[2]);;
//  e = 1e-6; //e is the tolerance
 // lambda = 0;
  
  
  

 /* MPI_Request req_A_send;
  MPI_Request req_A_recv;
*/
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  blockrows_size = mat_size / nprocs ;

  double A_blocks[mat_size*blockrows_size];
  generate_matrix(A_blocks,mat_size);
  /* Define new type to send datas */
 /* MPI_Datatype blocktype, blockselect;

  MPI_Type_contiguous(mat_size*blockrows_size, MPI_DOUBLE, &blocktype);
  MPI_Type_commit(&blocktype);

  MPI_Type_vector(blockrows_size, mat_size, mat_size, MPI_DOUBLE, &blockselect);
  MPI_Type_commit(&blockselect);
*/
 // if(rank == 0) {
    /*We initialize the matrix A and send rows of A to the other processes */
  //  A = (double *)malloc(sizeof(double) * mat_size * mat_size);

    
    /* Matrix generation*/
   /* for(i = 0; i < mat_size; i++){
       for(j = 0; j < mat_size; j++){
         // if (i == j){
             A[i*mat_size + j] = 1;
         // }*/
        /*  else if(i+1 == j && i < mat_size-1){
             A[i*mat_size + j] = mat_size;
          }
          else if(i == j+1 && j < mat_size-1){
             A[i*mat_size + j] = mat_size;
          }
          else{
             A[i*mat_size +j] = 0;
          }*/
     //  }    
   // }
   // Prvalues(mat_size, mat_size, A);
    
 
    /* Distribution of Block_rows of A */
    /* for (i = 0; i < nprocs; i++){
            int i_b_m = i * blockrows_size * mat_size;
            MPI_Isend(&A[i_b_m], 1, blockselect, i, 3, MPI_COMM_WORLD, &req_A_send);
            
     }*/
     
     /* double norm = compute_norm(2,x);
      printf("norm %.3f\t \n",norm);         // Testing the normalization of the vectors
      norm_vect(2,x);
      Prvalues(1,2,x);*/       
 // }

  

  /* After receiving rows of A, we start the power method and use Allgather to send y to all processes */
 /* A_blockrows = (double *)malloc(mat_size * blockrows_size * sizeof(double));
  MPI_Irecv(A_blockrows, 1, blocktype, 0, 3, MPI_COMM_WORLD, &req_A_recv);

  MPI_Wait(&req_A_recv, MPI_STATUS_IGNORE);
  */
   
   /*for(j = 0; j<mat_size; j++){
      for(i = 0; i<blockrows_size; i++){
     
         if (j-rank*blockrows_size>i){
             A_blocks[j + i*mat_size] = 0; 
         } 
         else{
             A_blocks[j + i*mat_size] = rank*(blockrows_size) + i + 1;
         }
       }
   } */ 
 
    
  double lambda = power_method(A_blocks,M,mat_size);
  
  

  /* Printing Part*/
  sleep(rank);
  printf("rank %d \n", rank);
  //Prvalues(mat_size,1,x);
  //Prvalues(blockrows_size,1,y_blocks);
  //Prvalues(mat_size,1,y);
  if (rank == 0) printf("lambda = %.2f\t \n",lambda);
  //if (rank == 2) Prvalues(1, blockrows_size, y);
  Prvalues(mat_size, blockrows_size, A_blocks);

  /*When the algorithm terminates we print the dominant eigein value and the associated eigen vector in root processor */ 

  /*free(y_blocks);
  free(y);
  free(x);
  free(A_blockrows);*/


 /* MPI_Type_free(&blockselect);
  MPI_Type_free(&blocktype);
  */
  MPI_Finalize();


  return 0;

}



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

void mat_mult( int length, double A_rows[], double x[], double y[]){

    int  i, j, rank, nprocs;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);


    /* Broadcast the column vector from proc 0 to all processes */
    MPI_Bcast(x,length,MPI_DOUBLE,0,MPI_COMM_WORLD);

    double tmp[length/nprocs];      
    
    for(i = 0; i < length/nprocs; i++){
        for(j = 0; j < length ; j++){
           tmp[i] += A_rows[i*length + j] * x[j] ;
        }
    }

    MPI_Gather(tmp,length/nprocs,MPI_DOUBLE,y,length/nprocs,MPI_DOUBLE,0,MPI_COMM_WORLD);    
}

double power_method(double A_blocks[],int n_iterations,int mat_size){

    int rank, i, n;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    double lambda = 0;
   // double lambda0 = lambda + 2*e ;
   // int k = 0 ;

    /* Define starting vector */
    double x[mat_size];
    if(rank==0){
      for(i = 0; i < mat_size; i++){
         x[i] = 1;
      }
    }

    
     
   // while (abs(lambda-lambda0) >= e && k <= n_iterations ){
    for(n = 0; n < n_iterations; n++){

        
        //Prvalues(mat_size,1,x);

        if (rank == 0){
        norm_vect(mat_size,x);
        }

       
        //Prvalues(mat_size,1,x);
        double  tmp[mat_size];
        mat_mult(mat_size, A_blocks, x, tmp);

                

        if(rank ==0){
        //lambda0 = lambda ;
          for(i = 0; i < mat_size; i++){
              x[i] = tmp[i] ;             
             }
        
        }
        
        
    }

  
        lambda = compute_norm(mat_size, x) ;  
        return lambda ;
}



void generate_matrix(double matrix[], int mat_size){
  
     int i,j, rank, nprocs;

     MPI_Comm_rank(MPI_COMM_WORLD,&rank);
     MPI_Comm_size(MPI_COMM_WORLD,&nprocs);

     for(j=0; j< mat_size; j++){
        for(i=0; i<mat_size/nprocs;i++){
            if(j-rank*(mat_size/nprocs)>i){
                 matrix[i*mat_size + j] = 0;
            }
            else{
                 matrix[i*mat_size + j] = rank*(mat_size/nprocs)+ i +1;
            }
        }
     }



}
