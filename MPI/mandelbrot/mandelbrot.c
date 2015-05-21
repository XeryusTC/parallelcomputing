/* file: mandelbrot.c
 *(C) Arnold Meijster and Rob de Bruin
 *
 * A simple program for computing images of the Mandelbrot set.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

#define WIDTH  4096
#define HEIGHT 3072

#define MAXITER 3000

int rank, size;
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
        double scale, Image image)
{ /* routine that computes an image of the mandelbrot set */
    int w=image->width, h=image->height, **im=image->imdata;
    double a, b;
    double x, y, z;
    int i, j, k;

    for (i = 0; i < h; i++)
    {
        b = centerY + i*scale - ((h/2)*scale);
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
            im[i][j] = k;
        }
    }
}

int main(int argc, char **argv)
{
    Image mandelbrot;

    /* Initialize MPI */
    MPI_Init (&argc, &argv);
    MPI_Comm_rank (MPI_COMM_WORLD, &rank);
    MPI_Comm_size (MPI_COMM_WORLD, &size);

    mandelbrot = makeImage(WIDTH, HEIGHT);
    mandelbrotSet (-0.65, 0, 2.5/HEIGHT, mandelbrot);

    writePPM (mandelbrot, "mandelbrot.ppm");

    freeImage(mandelbrot);

    MPI_Finalize ();
    return 0;
}
