

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <string.h>

void Prvalues(int length, int heigth,  double matrix[length * heigth]);
void mat_mult(double A_rows[], double vect[], double result[], int length);
double compute_norm2(double vect[], int vect_length);
double powerMethod(double A_rows[], int n_times, int mat_size);



void generate_triangle(double A_rows[], int mat_size);
void generate_diagonal(double A_rows[], int mat_size);
void generate_ones(double A_rows[], int mat_size);
void generate_tridiag(double A_rows[], int mat_size);


int main(int argc, char* argv[]) {

  /* Set parameters */
  int mat_size = atoi(argv[1]);
  int n_times = atoi(argv[2]);


  /* Start parallel calculations */
  int rank, nprocs;
  MPI_Init(&argc, & argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);


  /* p has to be a divisor of N !! */


 
  int blockrows_size = mat_size / nprocs;

  /* Generate the matrix */
  double A_blocks[mat_size*blockrows_size];
  generate_tridiag(A_blocks,mat_size);
  
 
  /* Run the powerMethod algorithm */
  double start = MPI_Wtime();
  double lambda = powerMethod(A_blocks,n_times,mat_size);
  double stop = MPI_Wtime();



  /* Calculating times */
  double timediff = stop - start;
  double times[nprocs];
  double average = 0;
  MPI_Gather(&timediff,1,MPI_DOUBLE,times,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
  if(rank==0){
    for(int i=0; i<nprocs; i++)
      average += times[i];
    average /= nprocs;
  }


  /* Print the results */
  if(rank==0)
    printf("dimension\t%d\tdominant lambda\t%f\ttime\t%f\n",mat_size,lambda,average);


  /* End MPI */
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
    result += pow(vect[i],2);

  /* The norm of the vector is the square root of the sum */
  return pow(result,0.5);
}



void mat_mult(double A_rows[], double vect[], double result[], int mat_size){

  int rank, nprocs;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);


  /* Broadcast the vector defined on pr0 */
  MPI_Bcast(vect, mat_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);


  /* Set temporary result to zero */
  double y_i[mat_size/nprocs];
  for(int i=0; i<mat_size/nprocs; i++)
    y_i[i]=0;


  /* Perform the actual matrix multiplication in parallel */
  for(int i=0; i<mat_size/nprocs; i++)
  for(int j=0; j<mat_size  ; j++)
    y_i[i] += A_rows[j+mat_size*i] * vect[j];


  /* Gather the results from the different processors */
  MPI_Gather(y_i,mat_size/nprocs,MPI_DOUBLE,result,mat_size/nprocs,MPI_DOUBLE,0,MPI_COMM_WORLD);

}



double powerMethod(double A_rows[], int n_times, int mat_size) {

  int rank;
  double lambda=0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  

  /* Define starting vector */
  double x[mat_size];
  if(rank==0){
    
    for(int i=0; i<mat_size; i++)
      x[i]=1;
  }


  /* The actual iteration process */
  for(int n=0; n<n_times; n++){


    /* Normalize the vector*/
    if(rank==0){
      double norm = compute_norm2(x,mat_size);
      for(int i=0; i<mat_size; i++)
	x[i] = x[i]/norm;
    }


    /* Calculate the product of the matrix and the vector */
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
