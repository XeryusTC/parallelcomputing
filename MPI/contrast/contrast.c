/* file: contrast.c
 *(C) Arnold Meijster and Rob de Bruin
 *
 * A simple program for contrast stretching PPM images.
 */

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>

int rank, size;

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
	im->imdata[0] = safeMalloc(h*w*sizeof(int));
	for (row = 1; row < h; row++)
	{
		im->imdata[row] = im->imdata[0] + row * w;
	}
	return im;
}

static void freeImage(Image im)
{ /* routine for deallocating memory occupied by an image */
	free(im->imdata[0]);
	free(im->imdata);
	free(im);
}

static Image readPGM(char *filename)
{ /* routine that returns an image that is read from a PGM file */
	int c, w, h, maxval, row, col;
	FILE *f;
	Image im;
	unsigned char *scanline;

	if ((f = fopen(filename, "rb")) == NULL)
	{
		error("Opening of PGM file failed\n");
	}
	/* parse header of image file (should be P5) */
	if ((fgetc(f) != 'P') || (fgetc(f) != '5') || (fgetc(f) != '\n'))
	{
		error("File is not a valid PGM file\n");
	}
	/* skip commentlines in file (if any) */
	while ((c=fgetc(f)) == '#')
	{
		while ((c=fgetc(f)) != '\n');
	}
	ungetc(c, f);
	/* read width, height of image */
	col = fscanf (f, "%d %d\n", &w, &h);
	/* read maximum greyvalue (dummy) */
	col = fscanf (f, "%d\n", &maxval);
	if (maxval > 255)
	{
		error ("Sorry, readPGM() supports 8 bits PGM files only.\n");
	}
	/* allocate memory for image */
	im = makeImage(w, h);
	/* read image data */
	scanline = malloc(w*sizeof(unsigned char));
	for (row = 0; row < h; row++)
	{
		col = fread(scanline, 1, w, f);
		for (col = 0; col < w; col++)
		{
			im->imdata[row][col] = scanline[col];
		}
	}
	free(scanline);
	fclose(f);
	return im;
}

static void writePGM(Image im, char *filename)
{ /* routine that writes an image to a PGM file.
   * This routine is only handy for debugging purposes, in case
   * you want to save images (for example template images).
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

static Image readImage (char *filename)
{
	Image im;
	im = readPGM(filename);
	return im;
}

static void writeImage (Image im, char *filename)
{
	writePGM(im, filename);
}

static void contrastStretch(int *data, int length, int low, int high, int min,
		int max)
{
	float scale = (float)(high - low) / (max - min);
	int i;

	for (i=0; i<length; ++i)
		data[i] = scale * data[i];
}

static void calculateMinMax(int *data, int length, int *minmax)
{
	int i;
	minmax[0] = 255;
	minmax[1] = 0;
	for (i=0; i<length; ++i)
	{
		minmax[0] = (data[i] < minmax[0] ? data[i] : minmax[0]);
		minmax[1] = (data[i] > minmax[1] ? data[i] : minmax[1]);
	}
}

static void masterProcess(char *infile, char *outfile)
{
	int portion, *sendcount, *offset, i, min, max, minmax[2], *data;
	Image im;

	im = readImage(infile);
	portion = ceil(im->height * im->width / (float)size);

	/* work distribution */
	printf("Calculating work distribution\n");
	sendcount = safeMalloc(sizeof(int)*size);
	offset    = safeMalloc(sizeof(int)*size);
	data      = safeMalloc(sizeof(int)*portion);
	for (i=0; i<size; ++i)
	{
		if (i == (size-1))
			sendcount[i] = (im->height * im->width) - ((size - 1) * portion);
		else
			sendcount[i] = portion;
		offset[i] = i * portion;
		/* Send size to slave */
		if (i > 0)
			MPI_Send(&sendcount[i], 1, MPI_INT, i, 0, MPI_COMM_WORLD);
	}
	/* distribute work */
	printf("Distributing work\n");
	MPI_Scatterv(im->imdata[0], sendcount, offset, MPI_INT, data, portion,
			MPI_INT, 0, MPI_COMM_WORLD);

	/* MPI_Allreduce doesn't support different length data sizes, so we wrote
	 * our own instead */
	calculateMinMax(data, portion, minmax);
	MPI_Allreduce(&minmax[0], &min, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
	MPI_Allreduce(&minmax[1], &max, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

	printf("Stretching contrast\n");
	contrastStretch(data, portion, 0, 255, min, max);
	/* Gather all the data and write the file */
	MPI_Gatherv(data, portion, MPI_INT, im->imdata[0], sendcount,
			offset, MPI_INT, 0, MPI_COMM_WORLD);
	printf("Writing file\n");
	writeImage(im, outfile);

	freeImage(im);
	free(sendcount);
	free(offset);
	free(data);
}

static void slaveProcess()
{
	int length, *data, minmax[2], min, max;
	MPI_Status status;

	MPI_Recv(&length, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
	data = safeMalloc(length * sizeof(int));
	MPI_Scatterv(NULL, NULL, NULL, MPI_INT, data, length, MPI_INT, 0,
			MPI_COMM_WORLD);
	calculateMinMax(data, length, minmax);
	MPI_Allreduce(&minmax[0], &min, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
	MPI_Allreduce(&minmax[1], &max, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
	contrastStretch(data, length, 0, 255, min, max);
	MPI_Gatherv(data, length, MPI_INT, NULL, NULL, NULL, MPI_INT, 0,
			MPI_COMM_WORLD);
}

int main(int argc, char **argv)
{
	double t, tmax;
	/* Initialize MPI */
	MPI_Init (&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	if (argc != 3)
	{
		printf ("Usage: %s input.pgm output.pgm\n", argv[0]);
		exit (EXIT_FAILURE);
	}

	t = MPI_Wtime();
	if (rank == 0)
		masterProcess(argv[1], argv[2]);
	else
		slaveProcess();

	t = MPI_Wtime() - t;
	printf("%2d time: %f\n", rank, t);
	MPI_Reduce(&t, &tmax, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	if (rank == 0)
		printf("Largest time: %f\n", tmax);

	MPI_Finalize ();
	return EXIT_SUCCESS;
}
