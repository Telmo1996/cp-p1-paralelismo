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

int main(int argc, char *argv[])
{
	MPI_Init(&argc, &argv);

    int i, done = 0, n;
    double PI25DT = 3.141592653589793238462643;
    double pi, h, sum, x;
	int numProcs, rank;
	int N[1];
	double Pi[1];

	MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	
	if(rank == 0){

		while (!done)
		{
			printf("Enter the number of intervals: (0 quits) \n");
			scanf("%d",&n);

			N[0] = n;
			for(i=1; i<numProcs; i++){
				MPI_Send(N, 1, MPI_INT, i, 99, MPI_COMM_WORLD);
			}

			if (n == 0) break;

			pi = calculaPi(0, numProcs, n);

			for(i=1; i<numProcs; i++){
				MPI_Recv(Pi, 1, MPI_DOUBLE, MPI_ANY_SOURCE, 99, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				pi += Pi[0];
			}

			printf("pi is approximately %.16f, Error is %.16f\n", pi, fabs(pi - PI25DT));
		}

	}else{
		while(!done){
			MPI_Recv(N, 1, MPI_INT, 0, 99, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			if (N[0] == 0) break;

			pi = calculaPi(rank, numProcs, N[0]);
			
			Pi[0] = pi;
			MPI_Send(Pi, 1, MPI_DOUBLE, 0, 99, MPI_COMM_WORLD);
		}

	}


	MPI_Finalize();
}
