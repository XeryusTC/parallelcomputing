/* file: mandelbrot.c
 *(C) Arnold Meijster and Rob de Bruin
 *
 * A simple program for computing images of the Mandelbrot set.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <string.h>

#define WIDTH 4096
#define HEIGHT 3072

#define MAXITER 3000

int rank, size, interval, start;
unsigned char colour[MAXITER][3];  /* colour look-up table */

/* Data type for storing 2D greyscale image */
typedef struct imagestruct
{
    int width, height;
    int **imdata;
} *Image;

static void error(char *errmsg)
{ /* print error message an abort program */
    fprintf (stderr, "%s\n", errmsg);
    exit(EXIT_FAILURE);
}

static void *safeMalloc(int n)
{ /* wrapper function for malloc with error checking */
    void *ptr = malloc(n);
    if (ptr == NULL)
    {
        error("Error: memory allocation failed.\n");
    }
    return ptr;
}

static Image makeImage(int w, int h)
{ /* routine for constructing (memory allocation) of images */
    Image im;
    int row;
    im = malloc(sizeof(struct imagestruct));
    im->width  = w;
    im->height = h;
    im->imdata = safeMalloc(h*sizeof(int *));
    for (row = 0; row < h; row++)
    {
        im->imdata[row] = safeMalloc(w*sizeof(int));
    }
    return im;
}

static void freeImage(Image im)
{ /* routine for deallocating memory occupied by an image */
    int row;
    for (row = 0; row < im->height; row++)
    {
        free(im->imdata[row]);
    }
    free(im->imdata);
    free(im);
}

static void writePPM(Image im, char *filename)
{ /* routine that writes an image to a PPM file. */
    int row, col;
    unsigned char *scanline;
    int i, h, idx;

    FILE *f = fopen(filename, "wb");
    if (f == NULL)
    {
        error("Opening of file failed\n");
    }
    /* write header of image file (P5) */
    fprintf(f, "P6\n");
    fprintf(f, "%d %d\n255\n", im->width, im->height);
    /* write image data */
    for (i=0; i<MAXITER; i++)
    {
        double angle = 3.14159265*i/MAXITER;
        colour[i][0] = (unsigned char)(128 + 127*sin(100*angle));
        colour[i][1] = (unsigned char)(128 + 127*sin(50*angle));
        colour[i][2] = (unsigned char)(128 + 127*sin(10*angle));
    }

    scanline = malloc(im->width*3*sizeof(unsigned char));
    for (row = 0; row < im->height; row++)
    {
        idx = 0;
        for (col = 0; col < im->width; col++)
        {
            h = im->imdata[row][col];
            scanline[idx++] = colour[h][0];
            scanline[idx++] = colour[h][1];
            scanline[idx++] = colour[h][2];
        }
        fwrite(scanline, 3, im->width, f);
    }
    free (scanline);
    fclose(f);
}

static void mandelbrotSet(double centerX, double centerY,
        double scale, Image image, int row)
{ /* routine that computes an image of the mandelbrot set */
    int w=image->width, h=image->height, **im=image->imdata;
    double a, b;
    double x, y, z;
    int j, k;

    b = centerY + row*scale - ((h/2)*scale);
    for (j = 0; j < w; j++)
    {
        a = centerX + j*scale - ((w/2)*scale);
        x = a;
        y = b;
        k = 0;
        while ((x*x + y*y <= 100) && (k < MAXITER))
        {
            z = x;
            x = x*x - y*y + a;
            y = 2*z*y + b;
            k++;
        };
        im[row][j] = k;
    }
}

int main(int argc, char **argv)
{
    Image mandelbrot;
    int i, running, data[WIDTH];
    MPI_Status status;
    double t;

    /* Initialize MPI */
    MPI_Init (&argc, &argv);
    MPI_Comm_rank (MPI_COMM_WORLD, &rank);
    MPI_Comm_size (MPI_COMM_WORLD, &size);

    mandelbrot = makeImage(WIDTH, HEIGHT);

    interval = mandelbrot->height / size;
    start = rank * interval;
    running = size - 1;


    /* reconstruct image */
    t = MPI_Wtime();
    if (rank == 0) {
        i = 0;
        while (running > 0) {
            MPI_Recv(data, mandelbrot->width, MPI_INT, MPI_ANY_SOURCE,
                    MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            /* slave is returning calculated data, add it to the image */
            if (status.MPI_TAG < mandelbrot->height) {
                memcpy(mandelbrot->imdata[status.MPI_TAG], data,
                        sizeof(int) * mandelbrot->width);
            }

            /* we are done doing calculations, send termination signal */
            if (i >= mandelbrot->height) {
                MPI_Send(&i, 1, MPI_INT, status.MPI_SOURCE, 0, MPI_COMM_WORLD);
                running--;
                printf("Slave %d shut down\n", status.MPI_SOURCE);
            /* send new data to calculate otherwise */
            } else {
                MPI_Send(&i, 1, MPI_INT, status.MPI_SOURCE, 1, MPI_COMM_WORLD);
                i++;
            }

        }

        writePPM (mandelbrot, "mandelbrot.ppm");
    } else {
        /* register with the master */
        MPI_Send(mandelbrot->imdata[0], mandelbrot->width, MPI_INT, 0,
                mandelbrot->height, MPI_COMM_WORLD);
        while (1) {
            MPI_Recv(&i, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            /* check if we received a kill signal */
            if (status.MPI_TAG == 0)
                break;
            mandelbrotSet (-0.65, 0, 2.5/HEIGHT, mandelbrot, i);
            MPI_Send(mandelbrot->imdata[i], mandelbrot->width, MPI_INT, 0, i,
                    MPI_COMM_WORLD);
        }
    }
    printf("process %d: %f\n", rank, MPI_Wtime()-t);

    freeImage(mandelbrot);

    MPI_Finalize ();
    return 0;
}
