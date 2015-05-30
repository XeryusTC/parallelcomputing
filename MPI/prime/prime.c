#include <stdio.h>
#include <math.h>
#include <mpi.h>

#define FALSE 0
#define TRUE  1

static int isPrime(unsigned int p)
{
	int i, root;
	if (p == 1)
		return FALSE;
	if (p == 2)
		return TRUE;
	if (p%2 == 0)
		return FALSE;

	root = (int)(1 + sqrt(p));
	for (i = 3; (i < root) && (p%i != 0); i += 2);
	return (i < root ? FALSE : TRUE);
}

int main (int argc, char **argv)
{
	unsigned int i, a, b, cnt=0, interval, message[2], tot;
	int size, rank;
	MPI_Status info;
	double t;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	/* MASTER */
	if (rank == 0)
	{
		printf ("Enter two integer number a, b such that 1<=a<=b: ");
		fflush(stdout);
		scanf ("%u %u", &a, &b);
		interval = (b - a) / size;

		message[0] = a;
		message[1] = a + interval;
		for (i=1; i<size; ++i)
		{
			message[0] += interval + 1;
			message[1] += interval + 1;
			if (i == size - 1)
				message[1] = b;
			MPI_Send(message, 2, MPI_INT, i, 0, MPI_COMM_WORLD);
		}
		b = a + interval;
	}
	/* SLAVE */
	else
	{
		MPI_Recv(message, 2, MPI_INT, 0, 0, MPI_COMM_WORLD, &info);
		a = message[0];
		b = message[1];
	}
	printf("%d: interval %d - %d\n", rank, a, b);

	/* Calculation code is shared */
	t = MPI_Wtime();
	if (a <= 2)
	{
		cnt = 1;
		a = 3;
	}
	if (a % 2 == 0)
	{
		a++;
	}

	for (i = a; i <= b; i += 2)
		if (isPrime(i))
			cnt++;
	printf("Time spend in %2d: %f\n", rank, MPI_Wtime() - t);

	MPI_Reduce(&cnt, &tot, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	if (rank == 0)
		printf ("\n#primes=%u\n", tot);

	MPI_Finalize();
	return 0;
}
