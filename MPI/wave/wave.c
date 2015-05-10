/* compilation: icc -Wall -mcmodel=large -O3 wave.c -lm  */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <time.h>
#include <math.h>

typedef float real;
typedef unsigned char byte;

int  N, NFRAMES;
int  nsrc, *src, *ampl;
real gridspacing, timespacing, speed;
int  iter;
real ***u;

#define ABS(a) ((a)<0 ? (-(a)) : (a))

static void initialize(real dx, real dt, real v, int n);
static void boundary(void);
static void solveWave(void);
static void stretchContrast(void);
static void saveFrames(int bw);
static void parseIntOpt(int argc, char **argv, char *str, int *val) ;
static void parseRealOpt(int argc, char **argv, char *str, real *val) ;


static void initialize(real dx, real dt, real v, int n)
{
  int i, j, k;
  real lambda;

  timespacing = dt;
  gridspacing = dx;
  speed = v;

  lambda = speed*timespacing/gridspacing;
  if (lambda > 0.5*sqrt(2))
  {
    printf ("Error: Convergence criterion is violated.\n");
    printf ("speed*timespacing/gridspacing=");
    printf ("%lf*%lf/%lf=%f>0.5*sqrt(2)\n", 
            speed, timespacing, gridspacing,lambda);
    timespacing=(real)(gridspacing*sqrt(2)/(2*speed));
    printf ("Timestep changed into: timespacing=%lf\n", timespacing);
  }

  /* draw n random locations of wave sources */
  srand(time(NULL));  /* initialize random generator with time */
  nsrc = n;
  src = malloc(2*nsrc*sizeof(int));
  ampl = malloc(nsrc*sizeof(int));
  for (i=0; i<nsrc; i++)
  {
    src[2*i] = random()%N;
    src[2*i+1] = random()%N;
    ampl[i] = 1;//1+random()%5;
  }

  /* allocate memory for u */
  u = malloc(NFRAMES*sizeof(real **));
  for (k=0; k<NFRAMES; k++)
  {
    u[k] = malloc(N*sizeof(real *));
    u[k][0] = malloc(N*N*sizeof(real *));
    for (i=1; i<N; i++)
    {
      u[k][i] = u[k][i-1] + N;
    }
  }

  /* initialize first two time steps */
  for (i=0; i<N; i++)
  {
    for (j=0; j<N; j++)
    {
      u[0][i][j] = 0;
    }
  }
  for (i=0; i<N; i++)
  {
    for (j=0; j<N; j++)
    {
      u[1][i][j] = 0;
    }
  }
}

static void boundary(void)
{
  int i, j;
  real t = iter*timespacing;

  for (i=0; i<N; i++)
  {
    for (j=0; j<N; j++)
    {
      u[iter][i][j] = 0;
    }
  }
  for (i=0; i<nsrc; i++)
  {
    u[iter][src[2*i]][src[2*i+1]] = (real)(ampl[i]*sin(t));
  }
}

static void solveWave(void)
{
  real sqlambda;
  int i, j;
  sqlambda = speed*timespacing/gridspacing;
  sqlambda = sqlambda*sqlambda;

  for (iter=2; iter<NFRAMES; iter++)
  {
    boundary();
    for (i=1; i<N-1; i++)
    {
      for (j=1; j<N-1; j++)
      {
	    u[iter][i][j] += 
	       sqlambda*(u[iter-1][i+1][j] + u[iter-1][i-1][j] + 
		             u[iter-1][i][j+1] + u[iter-1][i][j-1])
	       + (2-4*sqlambda)*u[iter-1][i][j]
	       - u[iter-2][i][j];
      }
    }
  }
}

static void stretchContrast(void)
{
  int i, j, frame;
  real min, max, scale = 255.0;
  
  min =  9999;
  max = -9999;
  for (frame=2; frame<NFRAMES; frame++)
  {
    for (i=0; i<N; i++)
    {
      for (j=0; j<N; j++)
      {
        if (u[frame][i][j] < min)
        {
          min = u[frame][i][j];
        }
        if (u[frame][i][j] > max)
        {
          max = u[frame][i][j];
        }
      }
    }
  }
  if (max>min)
  {
    scale = scale/(max-min);
  }
  for (frame=0; frame<NFRAMES; frame++)
  {
    for (i=0; i<N; i++)
    {
      for (j=0; j<N; j++)
      {
        u[frame][i][j] = scale*(u[frame][i][j] - min);
      }
    }
  }
}

static void saveFrames(int bw)
{
  FILE *ppmfile;
  char filename[32];
  int i, j, k, frame, idx;
  byte lut[256][3];  /* colour lookup table */
  byte *rgb;  /* framebuffer */

  if (bw)
  {
    /* grey values */
    for (i=0; i<256; i++)
    {
      lut[i][0] = lut[i][1] = lut[i][2] = (byte)i;
    }
  } else
  {
    /* color values */
    for (i=0; i<256; i++)
    {
      lut[i][0] = (byte)i;
      lut[i][1] = (byte)(127+2*(i<64 ? i : (i<192 ? 128-i : i-255)));
      lut[i][2] = (byte)(255-i);
    }
  }

  rgb = malloc(3*N*N*sizeof(byte));
  for (frame=0; frame<NFRAMES; frame++)
  {
    k = 0;
    for (i=0; i<N; i++)
    {
      for (j=0; j<N; j++)
      {
        idx = (int)u[frame][i][j];
        rgb[k++] = (byte)lut[idx][0];  /* red   */
        rgb[k++] = (byte)lut[idx][1];  /* green */
        rgb[k++] = (byte)lut[idx][2];  /* blue  */
      }
    }

    /* show color map (comment this out if you dont like it) */
    if (N >=255)
    {
      int i0 = N/2 - 128;
      for (i=0; i<256; i++)
      {
        k = 3*((i+i0)*N + 10);
        for (j=0; j<16; j++)
        {
          rgb[k++] = lut[i][0];
          rgb[k++] = lut[i][1];
          rgb[k++] = lut[i][2];
        }
      }
    }
    sprintf (filename, "frame%04d.ppm", frame);
    ppmfile = fopen(filename, "wb");
    fprintf (ppmfile, "P6\n");
    fprintf (ppmfile, "%d %d\n", N, N);
    fprintf (ppmfile, "255\n");
    fwrite (rgb, sizeof(byte), 3*N*N, ppmfile);
    fclose(ppmfile);
  }
  free(rgb);
}

void parseIntOpt(int argc, char **argv, char *str, int *val)
{
  int i, found=0;
  for (i=1; i<argc; i++) {
    if (strcmp(argv[i], str) == 0)
    {
      found++;
      if (found > 1)
      {
        printf ("Error: doubly defined option %s\n", str);
        exit(EXIT_FAILURE);
      }
      if (i+1 < argc)
      {
        *val = atoi(argv[i+1]);
      } else
      {
        printf ("Error: missing numeric  value after %s\n", str);
        exit(EXIT_FAILURE);
      }
    }
  }
}

void parseRealOpt(int argc, char **argv, char *str, real *val)
{
  int i, found=0;
  for (i=1; i<argc; i++)
  {
    if (strcmp(argv[i], str) == 0)
    {
      found++;
      if (found > 1)
      {
        printf ("Error: doubly defined option %s\n", str);
        exit(EXIT_FAILURE);
      }
      if (i+1 < argc)
      {
        *val = (real)atof(argv[i+1]);
      } else
      {
        printf ("Error: missing numeric  value after %s\n", str);
        exit(EXIT_FAILURE);
      }
    }
  }
}

int main (int argc, char **argv)
{
  int n;    /* number of sources                          */
  int bw;   /* colour(=0) or greyscale(=1) frames         */
  real dt;  /* time spacing (delta time)                  */
  real dx;  /* grid spacing (distance between grid cells) */
  real v;   /* velocity of waves                          */

  struct timeval start, end;
  double fstart, fend;

  /* default values */
  n = 10;
  bw = 0;
  NFRAMES = 100;   /* number of images/frames   */
  N  = 300;        /* grid cells (width,height) */
  dt = 0.1;
  dx = 0.1;
  v  = 0.5;

  /* any default setting changed by user ? */
  parseIntOpt(argc, argv, "-f", &NFRAMES);
  parseIntOpt(argc, argv, "-src", &n);
  parseIntOpt(argc, argv, "-bw", &bw);
  parseIntOpt(argc, argv, "-n", &N);
  parseRealOpt(argc, argv, "-t", &dt);
  parseRealOpt(argc, argv, "-g", &dx);
  parseRealOpt(argc, argv, "-s", &v);

  /* initialize system */
  initialize(dx, dt, v, n);

  /* start timer */
  printf ("Computing waves\n");
  gettimeofday (&start, NULL);

  /* solve wave equation */
  solveWave();

  /* stop timer */
  gettimeofday (&end, NULL);
  fstart = (start.tv_sec * 1000000.0 + start.tv_usec) / 1000000.0;
  fend = (end.tv_sec * 1000000.0 + end.tv_usec) / 1000000.0;
  printf ("wallclock: %lf seconds (ca. %5.2lf Gflop/s)\n", 
          fend-fstart, 
          (9.0*N*N*NFRAMES/(fend-fstart))/(1024*1024*1024));

  /* save images */
  printf ("Saving frames\n");
  stretchContrast();
  saveFrames(bw);

  printf ("Done\n");
  return EXIT_SUCCESS;
}
