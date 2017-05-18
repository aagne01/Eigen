CC      = mpicc

CCFLAGS = -O3 -Wall -g

LIBS    = -lmpi -lm

MPI_eigen: MPI_eigen.c
	$(CC) $(CCFLAGS) -o MPI_eigen MPI_eigen.c $(LIBS)

test_bcast: test_bcast.c
	$(CC) $(CCFLAGS) -o test_bcast test_bcast.c $(LIBS)

clean:
	$(RM) MPI_eigen
	$(RM) test_bcast
	$(RM) core.*
