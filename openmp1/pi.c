#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <sys/time.h>

#define f(x) (4.0/(1.0+(x)*(x)))

static double timer(void)
{
  struct timeval tm;
  gettimeofday (&tm, NULL);
  return tm.tv_sec + tm.tv_usec/1000000.0;
}

int main(int argc, char **argv)
{
  int i, N;
  double dx, sum, err;
  double clock;

  while (1)
  {
    /* Step (1): get a value for N */
    do {
      printf ("Enter number of approximation intervals:(0 to exit)\n");
      scanf ("%d", &N);
    } while (N<0);
    printf ("N=%d\n", N);

    /* Step (2): check for exit condition. */
    if (N == 0)
    {
      printf("Done\n");
      break;
    }

    /* start timer */
    clock = timer();

    /* Step (3): numerical integration */
    dx = 1.0/N;
    sum = 0.0;
    #pragma omp parallel for reduction(+:sum)
    for (i=0; i < N; i++)
    {
      sum = sum + f(dx*(i+0.5));
    }
    sum *= dx;

    /* stop timer */
    clock = timer() - clock;
    printf ("wallclock: %lf seconds\n", clock);

    /* Step (4): print the results */
    err = sum - 4*atan(1);
    printf ("%d: pi~%20.18lf   er~=%20.18lf\n", N, sum, err);
  }

  return EXIT_SUCCESS;
}
