#include <stdio.h>
#include <math.h>

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
  unsigned int i, a, b, cnt=0;

  printf ("Enter two integer number a, b such that 1<=a<=b: ");
  fflush(stdout);
  scanf ("%u %u", &a, &b);

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

  printf ("\n#primes=%u\n", cnt);
  return 0;
}
