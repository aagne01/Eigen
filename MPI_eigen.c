/*
* Parallel and Distributed Programming 
* Individual project
* author : Amadou Agne
* May 2017
*/

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <string.h>

void pr_values(int length, int heigth,  double matrix[length * heigth]);
void mat_mult(double A_rows[], double vect[], double result[], int length);
double compute_norm2(double vect[], int vect_length);
double power_method(double A_rows[], int n_times, int mat_size);

void generate_triangle(double A_rows[], int mat_size);
void generate_diagonal(double A_rows[], int mat_size);
void generate_ones(double A_rows[], int mat_size);
void generate_tridiag(double A_rows[], int mat_size);


int main(int argc, char* argv[]) {

  /* input parameters */
  int mat_size = atoi(argv[1]);
  int n_times = atoi(argv[2]);
  
  
  
  int rank, nprocs;
  MPI_Init(&argc, & argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
     
  /* nprocs has to be a divisor of mat_size ! */
  if (rank == 0){
     if(mat_size%nprocs != 0 || argc > 3 || argc < 2){
        printf("How to use: mpirun -np nprocs ./MPI_eigen mat_size ntimes \n");
        exit(0);
     }
   }
 
  int blockrows_size = mat_size / nprocs;

  /* Generate the matrix */
  double A_blocks[mat_size*blockrows_size];
  generate_ones(A_blocks,mat_size);
  
 
  /* Run the powerMethod algorithm */
  double start = MPI_Wtime();
  double lambda = power_method(A_blocks,n_times,mat_size);
  double stop = MPI_Wtime();



  /* Calculating times */
  double time_elapsed = stop - start;
  double times[nprocs];
  double mean = 0;
  MPI_Gather(&time_elapsed,1,MPI_DOUBLE,times,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
  if(rank==0){
    for(int i=0; i<nprocs; i++)
      mean += times[i];
    mean /= nprocs;
  }


  /* Print results */
  if(rank==0){
    printf("dominant lambda :  %.4f\t\n",lambda);
    printf("time of the power method algorithm :  %.4f\t\n",mean);
  }

 
  MPI_Finalize();

  return 0;
}


/* Function that generates a triangular matrix, used for verification purposes */
void generate_triangle(double A_rows[], int mat_size){

  int rank, nprocs;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  /* Define the matrix on the (different) processors */
  for(int j=0; j<mat_size  ; j++) {
  for(int i=0; i<mat_size/nprocs; i++) {

    if(j-rank*mat_size/nprocs>i)
      A_rows[j+mat_size*i] = 0;
    else
      A_rows[j+mat_size*i] = rank*(mat_size/nprocs)+i+1;

  }}

}

/*Function that generates a diagonal matrix, used for verification purposes*/

void generate_diagonal(double A_rows[], int mat_size){

  int rank, nprocs;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  /* Define the matrix on the (different) processors */
  for(int j=0; j<mat_size ; j++){
  for(int i=0; i<mat_size/nprocs; i++){
     if(j == i + rank*mat_size/nprocs){
       A_rows[j+mat_size*i] = i+1+rank*mat_size/nprocs;
     }

     else{
       A_rows[j+mat_size*i] = 0;
     }
      
  }
  }

}

/*Function that generates an all ones matrix, used for verification purposes*/

void generate_ones(double A_rows[], int mat_size){

  int rank, nprocs;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  /* Define the matrix on the (different) processors */
  for(int j=0; j<mat_size ; j++){
  for(int i=0; i<mat_size/nprocs; i++){
       A_rows[j+mat_size*i] = 1;
      
  }
  }

}


/*Function that generates a tri-diagonal matrix, used to measure performance*/

void generate_tridiag(double A_rows[], int mat_size){

  int rank, nprocs;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  /* Define the matrix on the (different) processors */
  for(int j=0; j<mat_size ; j++){
  for(int i=0; i<mat_size/nprocs; i++){
     if(j >1+ i + rank*mat_size/nprocs){
       A_rows[j+mat_size*i] = 0;
     }

     else if(j<i-1+rank*mat_size/nprocs){
       A_rows[j+mat_size*i] = 0;
     }

     else{
       A_rows[j+mat_size*i] = mat_size;
     }
      
  }
  }

}

double compute_norm2(double vect[], int vect_length){

  /* Add the square of each element */
  double result = 0;
  for(int i=0; i<vect_length;i++)
    result += vect[i]*vect[i];

  
  return sqrt(result);
}



void mat_mult(double A_rows[], double vect[], double result[], int mat_size){

  int rank, nprocs;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);


  /* Broadcast the vector x defined on proc 0 */
  MPI_Bcast(vect, mat_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);


  /* Set y_i  i=1 .. nprocs  to zero */
  double y_i[mat_size/nprocs];
  for(int i=0; i<mat_size/nprocs; i++)
    y_i[i]=0;


  /* Perform the actual matrix multiplication in parallel */
  for(int i=0; i<mat_size/nprocs; i++)
  for(int j=0; j<mat_size  ; j++)
    y_i[i] += A_rows[j+mat_size*i] * vect[j];


  /* Gather the results from the different processors and put it back to proc 0 */
  MPI_Gather(y_i,mat_size/nprocs,MPI_DOUBLE,result,mat_size/nprocs,MPI_DOUBLE,0,MPI_COMM_WORLD);

}



double power_method(double A_rows[], int n_times, int mat_size) {

  int rank;
  double lambda=0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  

  /* Defining and initializing vector x */
  double x[mat_size];
  if(rank==0){
    
    for(int i=0; i<mat_size; i++)
      x[i]=1;
  }


  /* The main loop of the algorithim */
  for(int n=0; n<n_times; n++){


    /* Normalizing vector x*/
    if(rank==0){
      double norm = compute_norm2(x,mat_size);
      for(int i=0; i<mat_size; i++)
	x[i] = x[i]/norm;
    }


    /* Computing the multiplication between the matrix rows and x */
    double x_tmp[mat_size];
    mat_mult(A_rows,x,x_tmp,mat_size);
    if(rank==0){
      for(int i=0; i<mat_size; i++)
	x[i] = x_tmp[i];
      }
      lambda = compute_norm2(x,mat_size);
    }

  /* The result is the norm of the vector */
  
  return lambda;
}




/* Function to Print matrix and vectors*/
void pr_values(int length, int heigth,  double matrix[length * heigth]){   
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
