/* compilation: icc -Wall -mcmodel=large -O3 wave.c -lm  */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <time.h>
#include <math.h>
#include <mpi.h>

typedef float real;
typedef unsigned char byte;

int  N, NFRAMES;
int  nsrc, *src, *ampl;
real gridspacing, timespacing, speed;
int  iter;
real ***u;
int rank, size;

#define ABS(a) ((a)<0 ? (-(a)) : (a))

static void initialize(real dx, real dt, real v, int n);
static void boundary(void);
static void solveWave(int *sendcount, int *offset);
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
    if (lambda > 0.5*sqrt(2) && rank == 0)
    {
        printf ("Error: Convergence criterion is violated.\n");
        printf ("speed*timespacing/gridspacing=");
        printf ("%lf*%lf/%lf=%f>0.5*sqrt(2)\n",
                speed, timespacing, gridspacing,lambda);
        timespacing=(real)(gridspacing*sqrt(2)/(2*speed));
        printf ("Timestep changed into: timespacing=%lf\n", timespacing);
    }
    else if (lambda > 0.5*sqrt(2))
    {
        /* slaves can change the timespacing quietly */
        timespacing=(real)(gridspacing*sqrt(2)/(2*speed));
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

static void generateSources(int n)
{
    int i;
    srand(time(NULL));

    nsrc = n;
    src = malloc(2 * nsrc * sizeof(int));
    ampl = malloc(nsrc * sizeof(int));

    for (i=0; i<nsrc; i++)
    {
        src[2*i] = random() % N;
        src[2*i+1] = random() % N;
        ampl[i] = 1;
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

static void solveWave(int *sendcount, int *offset)
{
    real sqlambda, *data;
    int i, j, start, end, *sizes, *offsets;

    sqlambda = speed*timespacing/gridspacing;
    sqlambda = sqlambda*sqlambda;
    start = offset[rank];
    end = offset[rank] + sendcount[rank];
    sizes   = malloc(sizeof(int) * size);
    offsets = malloc(sizeof(int) * size);
    for (i=0; i<size; ++i)
    {
        sizes[i]   = sendcount[i]*N;
        offsets[i] = offset[i]*N;
        if (rank == 0) {
            printf("Process %d sizes: %d\t%d\n", i, sizes[i], offsets[i]);
        }
    }
    data = malloc(sizeof(real)*sizes[rank]);

    for (iter=2; iter<NFRAMES; iter++)
    {
        boundary();
        for (i=start+1; i<end-1; i++)
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
        memcpy(data, (u[iter][0] + offsets[rank]), sizes[rank]*sizeof(real));
        MPI_Allgatherv(data, sizes[rank], MPI_FLOAT,
                u[iter][0], sizes, offsets, MPI_FLOAT, MPI_COMM_WORLD);
    }
    free(data);
    free(sizes);
    free(offsets);
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

static void masterProcess(int argc, char **argv)
{
    double t;
    int i, *sendcount, *offset, chunksize;
    int n = 10, bw = 0, intopts[4];
    real dt = 0.1, dx = 0.1, v = 0.5, realopts[3];
    NFRAMES = 100;
    N = 300;

    /* any default setting changed by user ? */
    parseIntOpt(argc, argv, "-f", &NFRAMES);
    parseIntOpt(argc, argv, "-src", &n);
    parseIntOpt(argc, argv, "-bw", &bw);
    parseIntOpt(argc, argv, "-n", &N);
    parseRealOpt(argc, argv, "-t", &dt);
    parseRealOpt(argc, argv, "-g", &dx);
    parseRealOpt(argc, argv, "-s", &v);

    generateSources(n);

    /* share the initial settings between everyone */
    intopts[0] = NFRAMES;
    intopts[1] = n;
    intopts[2] = bw;
    intopts[3] = N;
    realopts[0] = dt;
    realopts[1] = dx;
    realopts[2] = v;
    MPI_Bcast(intopts, 4, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(realopts, 3, MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Bcast(src, 2*nsrc, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(ampl, nsrc, MPI_INT, 0, MPI_COMM_WORLD);

    /* notify all the slaves of which part they have to do */
    chunksize = ceil(N / (float)size);
    sendcount = malloc(size * sizeof(int));
    offset    = malloc(size * sizeof(int));
    for (i=0; i<size; ++i)
    {
        if (i == (size-1))
            sendcount[i] = (N) - ((size - 1) * chunksize);
        else
            sendcount[i] = chunksize;
        offset[i] = i * chunksize;
    }
    MPI_Bcast(sendcount, size, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(offset, size, MPI_INT, 0, MPI_COMM_WORLD);
    for(i=0; i<size; ++i)
    {
        printf("Work for %2d: start %d\tsize %d\n", i, offset[i], sendcount[i]);
    }

    /* initialize system */
    initialize(dx, dt, v, n);

    printf("Solving wave equation...\n");
    t = MPI_Wtime();
    solveWave(sendcount, offset);

    /* save images */
    printf ("Saving frames\n");
    stretchContrast();
    saveFrames(bw);

    printf("%2d: elapsed time %f\n", rank, MPI_Wtime() - t);
}

static void slaveProcess()
{
    double t;
    int *sendcount, *offset;
    int n, bw, intopts[4];
    real dt, dx, v, realopts[3];
    MPI_Bcast(intopts, 4, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(realopts, 3, MPI_FLOAT, 0, MPI_COMM_WORLD);
    NFRAMES = intopts[0];
    n       = intopts[1];
    bw      = intopts[2];
    N       = intopts[3];
    dt      = realopts[0];
    dx      = realopts[1];
    v       = realopts[2];
    nsrc = n;
    src = malloc(2 * nsrc * sizeof(int));
    ampl = malloc(nsrc * sizeof(int));
    MPI_Bcast(src, 2*n, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(ampl, n, MPI_INT, 0, MPI_COMM_WORLD);
    /* receive which portion we have to calculate */
    sendcount = malloc(size *sizeof(int));
    offset = malloc(size *sizeof(int));
    MPI_Bcast(sendcount, size, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(offset, size, MPI_INT, 0, MPI_COMM_WORLD);

    /* initialize system */
    initialize(dx, dt, v, n);

    t = MPI_Wtime();
    solveWave(sendcount, offset);
    printf("%2d: elapsed time %f\n", rank, MPI_Wtime() - t);
}

int main (int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    printf("Process %2d reporting!\n", rank);

    if (rank == 0)
    {
        masterProcess(argc, argv);
        printf("Done\n");
    }
    else
    {
        slaveProcess();
    }
    MPI_Barrier(MPI_COMM_WORLD);
    return EXIT_SUCCESS;
}
