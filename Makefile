CC      = mpicc

CCFLAGS = -O3 -Wall -g

LIBS    = -lmpi -lm

MPI_eigen: MPI_eigen.c
	$(CC) $(CCFLAGS) -o MPI_eigen MPI_eigen.c $(LIBS)

clean:
	$(RM) MPI_eigen
