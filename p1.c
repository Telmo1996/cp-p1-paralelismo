#include <stdio.h>
#include <math.h>
#include <mpi.h>

double calculaPi(int numProc, int numProcs, int n){
	double h, sum, x;

	h   = 1.0 / (double) n;
	sum = 0.0;
	for (int i = numProc; i <= n; i += numProcs) {
		x = h * ((double)i - 0.5);
		sum += 4.0 / (1.0 + x*x);
	}
	return h * sum;
}

int MPI_BinomialColectiva(void *buffer, int count, MPI_Datatype datatype,
	int root, MPI_Comm comm){
	int i=1, rank, numprocs;
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	while (1){
		if(rank < pow(2,i-1)){
			if(rank+pow(2, i-1)>=numprocs){
				break;
			}
			MPI_Send(buffer, count, datatype, rank+pow(2, i-1), 0, comm);
			
		}else{
			if(rank-pow(2, i-1)<pow(2, i-1))
			MPI_Recv(buffer, count, datatype, rank-pow(2, i-1), 0, comm, MPI_STATUS_IGNORE);
		}
		i++;
	}	
	return 0;
}

int main(int argc, char *argv[])
{
	MPI_Init(&argc, &argv);

    int i, done = 0, n;
    double PI25DT = 3.141592653589793238462643;
    double pi, h, sum, x;
	int numProcs, rank;
	double pi_parcial;

	MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	
	if(rank == 0){

		while (!done)
		{
			printf("Enter the number of intervals: (0 quits) \n");
			scanf("%d",&n);
			
			MPI_BinomialColectiva(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

			if (n == 0) break;

			pi = calculaPi(0, numProcs, n);

			MPI_Reduce(&pi, &pi_parcial, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

			pi = pi_parcial;

			printf("pi is approximately %.16f, Error is %.16f\n", pi, fabs(pi - PI25DT));
		}

	}else{
		while(!done){

			MPI_BinomialColectiva(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
			if (n == 0) break;

			pi = calculaPi(rank, numProcs, n);
			
			MPI_Reduce(&pi, &pi_parcial, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		}

	}


	MPI_Finalize();
}
