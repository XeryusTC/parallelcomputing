#include <stdio.h>
#include <math.h>
#include <mpi.h>

#define FALSE 0
#define TRUE  1
#define BLOCK_SIZE pow(2,10)

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
	unsigned int i, j, a, b, cnt=0, interval, message[2];
	int size, rank;
	MPI_Status info;
	double t;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if (size > 0)
	{
		/* MASTER */
		if (rank == 0)
		{
			int tag[size], slaveid;
			for (i=0; i<size; ++i)
				tag[i] = 0;
			printf ("Enter two integer number a, b such that 1<=a<=b: ");
			fflush(stdout);
			scanf ("%u %u", &a, &b);
			/* The interval is smaller than the total when all slaves get an equal block */
			if (size * BLOCK_SIZE > b - a)
			{
				interval = (b - a) / size;
				/* Make sure every slave gets a block */
				for (i=1; i<size; ++i)
				{
					message[0] += interval + 1;
					message[1] += interval + 1;
					if (message[1] > b)
						message[1] = b;
					MPI_Send(message, 2, MPI_INT, i, tag[i], MPI_COMM_WORLD);
					tag[i]++;
				}
			}
			/* There are more blocks than slaves */
			else
			{
				message[0] = a;
				message[1] = a + BLOCK_SIZE;
				/* Send all the available blocks */
				while (message[1] < b)
				{
					slaveid = (cnt % (size - 1)) + 1;
					MPI_Send(message, 2, MPI_INT, slaveid, tag[slaveid], MPI_COMM_WORLD);
					message[0] += BLOCK_SIZE + 1;
					message[1] += BLOCK_SIZE + 1;
					if (message[1] > b)
						message[1] = b;
					cnt++;
					tag[slaveid] += 1;
				}
			}
			/* Send stop command */
			message[0] = message[1] = 0;
			for (i=1; i<size; ++i)
				MPI_Send(message, 2, MPI_INT, i, tag[i], MPI_COMM_WORLD);
		}
		/* SLAVE */
		else
		{
			j = 0;
			t = MPI_Wtime();
			while(TRUE)
			{
				MPI_Recv(message, 2, MPI_INT, 0, j, MPI_COMM_WORLD, &info);
				a = message[0];
				b = message[1];
				if (!a && !b)
					break;
				j++;

				/* Calculation code is shared */
				if (a <= 2)
				{
					cnt++;
					a = 3;
				}
				if (a % 2 == 0)
				{
					a++;
				}

				for (i = a; i <= b; i += 2)
					if (isPrime(i))
						cnt++;
			}
			printf("Time spend in %2d: %f\n", rank, MPI_Wtime() - t);
			MPI_Send(&cnt, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
		}
	}

	/* MASTER */
	if (rank == 0)
	{
		for (i=1; i<size; ++i)
		{
			MPI_Recv(message, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &info);
			cnt += message[0];
		}
		printf ("\n#primes=%u\n", cnt);
	}
	/* SLAVE */
	else
	{
	}

	MPI_Finalize();
	return 0;
}
