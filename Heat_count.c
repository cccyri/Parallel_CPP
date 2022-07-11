#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <mpi.h>

double pi = 3.1415926535897931159979635;

int main(int argc, char* argv[]) 
{
	int rank, size; //processes init

	double a, b, l;
	a = 0;
	b = 0;
	l = 1;

	int k = 1;
	int N = atoi(argv[1]) + 1; //correct dots init


	double h = l / (N - 1); //step
	double tau = 0.5 * h * h;
	double T = 500 * tau; //time

	MPI_Init(&argc, &argv); 
	MPI_Comm_size(MPI_COMM_WORLD, &size); //amt proc
	MPI_Comm_rank(MPI_COMM_WORLD, &rank); //rank proc


	//exact solution by formula
	if (rank == 0) 
	{
		double formula_sol = 0;
		printf("exact solution by formula:\n");
		for (int i = 0; i < 11; i++) 
		{
			formula_sol = 0;
			for (int j = 0; j < 100000; j++) 
			{
				formula_sol += 4 / pi * exp(-pi * pi * T * (2 * j + 1) * (2 * j + 1) / (l * l)) * sin(pi * 0.1 * i * (2 * j + 1) / l) / (2 * j + 1);
			}
			printf("u(%lf, %lf) = %lf\n", 0.1 * i, T, formula_sol);
		}
		printf("\n\n");

	}

	MPI_Barrier(MPI_COMM_WORLD);

	double start_time = MPI_Wtime(); //new time start

	double* Array_;
	double* Array_New;
	int MiniSize = N / size;

	if (rank == size - 1 && N % size != 0)
	{
		MiniSize += N % size;
	}
	
	int FullSize = MiniSize + 2;

	//memory for main elements and 2 places for boundary conditions
	Array_New = (double*)malloc((FullSize) * sizeof(double));
	Array_ = (double*)malloc((FullSize) * sizeof(double));
	int j;

	for (int i = 0; i < FullSize; i++)
	{
		Array_[i] = 1;
	}

	//first proc:
	if (rank == 0)
	{
		Array_[1] = a;
	}

	//last proc:
	if (rank == size - 1)
	{
		Array_[MiniSize] = b;
	}

	printf("rg: %d, h: %lf, T: %lf, tau: %10.3e T/tau: %d, N: %d, MiniSize: %d\n", rank, h, T, tau, (int)(T / tau), N, MiniSize);

	for (int i = 0; i < (int)(T / tau); i++) 
	{
		//geometrical paralleling v.3
		if (rank % 2 == 0)
		{
			if (rank > 0) 
			{
				MPI_Send(Array_ + 1, 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);
				MPI_Recv(Array_, 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			}
			if (rank < size - 1) 
			{
				MPI_Send(Array_ + MiniSize, 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
				MPI_Recv(Array_ + MiniSize + 1, 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			}
		}
		else 
		{
			if (rank > 0) 
			{
				MPI_Recv(Array_, 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				MPI_Send(Array_ + 1, 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);
			}
			if (rank < size - 1) 
			{
				MPI_Recv(Array_ + MiniSize + 1, 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				MPI_Send(Array_ + MiniSize, 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
			}
		}


		int left;
		if (rank == 0)
		{
			left = 2;
		}
		else left = 1;

		int right;
		if (rank == size - 1)
		{
			right = MiniSize;
		}
		else right = MiniSize + 1;


		for (int j = left; j < right; j++) 
		{
			Array_New[j] = Array_[j] + tau * k * (Array_[j - 1] - 2 * Array_[j] + Array_[j + 1]) / (h * h);
		}
		for (j = left; j < right; j++)
		{
			Array_[j] = Array_New[j];
		}
	}

	double x = rank * (N / size) * h;

	int left, right;

	left = 1;
	right = MiniSize + 1;

	for (int j = left; j < right; j++)
	{
		if (fabs((x * 10) - round(x * 10)) < 5 * h) 
		{
			printf("u(%lf, %lf) = %lf\n", x, T, Array_[j]);
		}
		x += h;
	}

	free(Array_);
	free(Array_New);

	MPI_Barrier(MPI_COMM_WORLD);

	double time = MPI_Wtime() - start_time;

	if (rank == 0)
	{
		printf("%lf\n", time);
	}

	MPI_Finalize();

	return 0;
}

