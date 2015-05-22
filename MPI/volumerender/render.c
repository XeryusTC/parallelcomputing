/* file: render.c
 *(C) Arnold Meijster and Rob de Bruin
 *
 * A simple orthogonal maximum intensity projection volume render.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

typedef unsigned char byte;

#define NFRAMES 360

int size, rank;

/* Data type for storing 2D greyscale image */
typedef struct imagestruct
{
    int width, height;
    int **imdata;
} *Image;

/* Data type for storing 3D greyscale volumes */
typedef struct volumestruct
{
    int width, height, depth;
    byte ***voldata;
    byte **voldata2d;
    byte *voldata3d;
} *Volume;

static void error(char *errmsg)
{ /* print error message an abort program */
    fprintf (stderr, "%s\n", errmsg);
    exit(EXIT_FAILURE);
}

static void *safeMalloc(int n)
{
    /* routine for safe memory allocation */
    void *ptr;
    ptr = malloc(n);
    if (ptr == NULL)
    {
        error("Error: malloc(%d) failed\n");
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
    im->imdata[0] = safeMalloc(h*w*sizeof(int));
    for (row = 1; row < h; row++)
    {
        im->imdata[row] = im->imdata[0] + row * w;
    }
    return im;
}

static void freeImage(Image im)
{ /* routine for deallocating memory occupied by an image */
    free(im->imdata);
    free(im);
}

static Volume makeVolume(int w, int h, int d)
{ /* routine for constructing (memory allocation) 3D volumes */
    Volume vol;
    int slice, row, shw;
    vol = malloc(sizeof(struct volumestruct));
    vol->width  = w;
    vol->height = h;
    vol->depth  = d;

    vol->voldata = safeMalloc(d*sizeof(byte**));
    vol->voldata2d = safeMalloc(d*h*sizeof(byte*));
    vol->voldata3d = safeMalloc(d*h*w*sizeof(byte));
    for (slice = 0; slice < d; slice++)
    {
        shw = slice * h * w;
        vol->voldata[slice] = vol->voldata2d + h * slice;
        for (row = 0; row < h; row++)
        {
            vol->voldata[slice][row] = vol->voldata3d + shw + row * w;
        }
    }
    return vol;
}

static void freeVolume(Volume vol)
{ /* routine for deallocating memory occupied by a volume */
    free(vol->voldata3d);
    free(vol->voldata2d);
    free(vol->voldata);
    free(vol);
}

static void writePGM(Image im, char *filename)
{ /* routine that writes an image to a PGM file.
*/
    int row, col;
    unsigned char *scanline;
    FILE *f = fopen(filename, "wb");
    if (f == NULL)
    {
        error("Opening of file failed\n");
    }
    /* write header of image file (P5) */
    fprintf(f, "P5\n");
    fprintf(f, "%d %d\n255\n", im->width, im->height);
    /* write image data */
    scanline = malloc(im->width*sizeof(unsigned char));
    for (row = 0; row < im->height; row++)
    {
        for (col = 0; col < im->width; col++)
        {
            scanline[col] = (unsigned char)(im->imdata[row][col] % 256);
        }
        fwrite(scanline, 1, im->width, f);
    }
    free (scanline);
    fclose(f);
}

static Volume readVolume (char *filename)
{ /* routine that reads in a volume data set from a VOX file. */
    Volume volume;
    byte ***vol;
    int i, j, k, width, height, depth;
    FILE *f;
    byte *scanline;

    f = fopen(filename, "rb");
    if (!f) error("Error opening volume file");

    k = fscanf (f, "D3\n%d %d %d\n255\n", &width, &height, &depth);

    volume = makeVolume(width, height, depth);
    vol = volume->voldata;
    scanline = safeMalloc(width*sizeof(byte));
    for (i = 0; i < depth; i++)
    {
        for (j = 0; j < height; j++)
        {
            k = fread(scanline, sizeof(byte), width, f);
            for (k = 0; k < width; k++)
            {
                vol[i][j][k] = scanline[k];
            }
        }
    }
    free(scanline);
    fclose (f);
    return volume;
}

static void rotateVolume (double rotx, double roty, double rotz,
        Volume volume, Volume rotvolume)
{
    /* rotate the volume around the x-axis with angle rotx, followed
     * by a rotation around the y-axis with angle roty, and finally
     * around the z-axis with angle rotz. The rotated volume is
     * returned in rotvolume.
     */
    int i, j, k, xi, yi, zi;
    int width  = volume->width;
    int height = volume->height;
    int depth  = volume->depth;
    byte ***vol = volume->voldata;
    byte ***rot = rotvolume->voldata;
    double x, y, z;
    double sinx, siny, sinz;
    double cosx, cosy, cosz;

    for (i=0; i<depth; i++)
    {
        for (j=0; j<height; j++)
        {
            for (k=0; k<width; k++)
            {
                rot[i][j][k] = 0;
            }
        }
    }

    sinx = sin(rotx);  siny = sin(roty); sinz = sin(rotz);
    cosx = cos(rotx);  cosy = cos(roty); cosz = cos(rotz);
    for (i=0; i<depth; i++)
    {
        for (j=0; j<height; j++)
        {
            for (k=0; k<width; k++)
            {
                xi = j - height/2;
                yi = k - width/2;
                zi = i - depth/2;

                /* rotation around x-axis */
                x = (double)xi;
                y = (double)(yi*cosx + zi*sinx);
                z = (double)(zi*cosx - yi*sinx);
                xi = (int)x;
                yi = (int)y;
                zi = (int)z;
                /* rotation around y-axis */
                x = (double)(xi*cosy + zi*siny);
                y = (double)(yi);
                z = (double)(zi*cosy - xi*siny);
                xi = (int)x;
                yi = (int)y;
                zi = (int)z;

                /* rotation around z-axis */
                x = (double)(xi*cosz + yi*sinz);
                y = (double)(yi*cosz - xi*sinz);
                z = (double)zi;

                xi = (int)(x + height/2);
                yi = (int)(y + width/2);
                zi = (int)(z + depth/2);
                if ((xi>=0) && (xi<height) &&
                        (yi>=0) && (yi<width) &&
                        (zi>=0) && (zi<depth))
                {
                    rot[zi][xi][yi] = vol[i][j][k];
                }
            }
        }
    }
    return;
}

static void contrastStretch (int low, int high, Image image)
{
    /* stretch the dynamic rabge of the image to the range [low..high] */
    int row, col, min, max;
    int width=image->width, height=image->height, **im=image->imdata;
    double scale;

    /* Determine minimum and maximum */
    min = max = im[0][0];
    for (row = 0; row < height; row++)
    {
        for (col = 0; col < width; col++)
        {
            min = (im[row][col] < min ? im[row][col] : min);
            max = (im[row][col] > max ? im[row][col] : max);
        }
    }

    /* Compute scale factor */
    scale = (double)(high-low) / (max-min);

    /* Stretch image */
    for (row = 0; row < height; row++)
    {
        for (col = 0; col < width; col++)
        {
            im[row][col] = (int)(scale*(im[row][col] - min));
        }
    }
}

static void orthoGraphicRenderer(Volume volume, Image image)
{
    /* Render image from volume (othographic maximum intensity
     * projection).
     */
    int i, j, k;
    int width=volume->width;
    int height=volume->height;
    int depth=volume->depth;
    int **im=image->imdata;
    byte ***vol=volume->voldata;

    for (i = 0; i < height; i++)
    {
        for (j = 0; j < width; j++)
        {
            im[i][j] = 0;
        }
    }

    for (i=0; i<depth; i++)
    {
        for (j=0; j<height; j++)
        {
            for (k=0; k<width; k++)
            {
                im[j][k] += vol[i][j][k];
            }
        }
    }

    contrastStretch(0, 255, image);
}

static void smoothImage(Image image, Image smooth)
{
    int width=image->width, height=image->height;
    int **im=image->imdata, **sm=smooth->imdata;
    int i, j, ii, jj, sum, cnt;

    for (i=0; i<height; i++)
    {
        for (j=0; j<width; j++)
        {
            cnt = 0;
            sum = 0;
            for (ii=i-1; ii<=i+1; ii++)
            {
                if ((ii>=0) && (ii<height))
                {
                    for (jj=j-1; jj<=j+1; jj++)
                    {
                        if ((jj >= 0) && (jj < width) && (im[ii][jj] != 0))
                        {
                            cnt++;
                            sum += im[ii][jj];
                        }
                    }
                }
            }
            sm[i][j] = (cnt == 0 ? 0 : sum / cnt);
        }
    }
}

void computeFrame2(int frame, double rotx, double roty, double rotz,
        Volume vol, Volume rot, Image image, Image smooth)
{
    rotateVolume(rotx, roty, rotz, vol, rot);
    orthoGraphicRenderer(rot, image);
    smoothImage(image, smooth);
}

int main (int argc, char **argv)
{
    int width, height, depth, v[3];
    Volume vol, rot;
    Image im, smooth;
    int frame, running, count;
    double rotx, roty, rotz, t;
    MPI_Status status;
    char fnm[256];

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    printf("Process %d reporting for duty!\n", rank);

    if (rank == 0)
    {
        if (argc!=2)
        {
            fprintf (stderr, "Usage: %s <volume.vox>\n", argv[0]);
            exit(EXIT_FAILURE);
        }

        vol = readVolume (argv[1]);
        width  = v[0] = vol->width;
        height = v[1] = vol->height;
        depth  = v[2] = vol->depth;
        im = makeImage(width, height);

        /* broadcast size to all processes so they can allocate memory */
        MPI_Bcast(v, 3, MPI_INT, 0, MPI_COMM_WORLD);
        /* send volume data to all processes */
        MPI_Bcast(vol->voldata3d, width*height*depth, MPI_BYTE, 0, MPI_COMM_WORLD);

        /* send frames to all the slaves */
        t = MPI_Wtime();
        running = size;
        count = 0;
        while (running > 1)
        {
            /* slave sends the frame which it is returning */
            MPI_Recv(&frame, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG,
                    MPI_COMM_WORLD, &status);
            /* retrieve the frame data from the slave */
            if (frame >= 0)
            {
                MPI_Recv(im->imdata[0], width*height, MPI_INT,
                        status.MPI_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD,
                        &status);
                /* write the frame to file */
                sprintf (fnm, "frame%04d.pgm", frame);
                writePGM(im, fnm);
            }

            if (count < NFRAMES)
            {
                /* send a new frame to the slave if any are available */
                MPI_Send(&count, 1, MPI_INT, status.MPI_SOURCE, 0,
                        MPI_COMM_WORLD);
                count++;
            }
            else
            {
                /* send termination signal */
                frame = -1;
                MPI_Send(&frame, 1, MPI_INT, status.MPI_SOURCE, 0,
                        MPI_COMM_WORLD);
                running--;
            }
        }
    }
    /* slave */
    else
    {
        /* retrieve the size so we can allocate memory */
        MPI_Bcast(v, 3, MPI_INT, 0, MPI_COMM_WORLD);
        width  = v[0];
        height = v[1];
        depth  = v[2];
        vol = makeVolume(width, height, depth);
        /* retrieve the volume data and store it in the volume */
        MPI_Bcast(vol->voldata3d, width*height*depth, MPI_BYTE, 0, MPI_COMM_WORLD);

        rot = makeVolume(width, height, depth);
        im = makeImage(width, height);
        smooth = makeImage(width, height);

        /* register with master */
        t = MPI_Wtime();
        frame = -1;
        MPI_Send(&frame, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
        while (1)
        {
            MPI_Recv(&frame, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD,
                    &status);
            /* shut down on term signal */
            if (frame == -1)
                break;
            /* compute frame and send it to the master */
            rotx = roty = rotz = 0;
            switch (3*frame/NFRAMES)
            {
                case 0: rotx = 6*M_PI*frame/NFRAMES; break;
                case 1: roty = 6*M_PI*frame/NFRAMES; break;
                case 2: rotz = 6*M_PI*frame/NFRAMES; break;
            }
            computeFrame2(frame, rotx, roty, rotz, vol, rot, im, smooth);
            MPI_Send(&frame, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
            MPI_Send(smooth->imdata[0], width*height, MPI_INT, 0, 0,
                    MPI_COMM_WORLD);
        }
        freeImage(smooth);
        freeVolume(rot);
    }
    printf("Process %d: %f sec\n", rank, MPI_Wtime()-t);
    freeVolume(vol);
    if (im != NULL);
        freeImage(im);

    MPI_Finalize();

    return 0;
}
