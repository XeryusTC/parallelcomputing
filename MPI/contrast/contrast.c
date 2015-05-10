/* file: contrast.c
 *(C) Arnold Meijster and Rob de Bruin
 *
 * A simple program for contrast stretching PPM images. 
 */

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

/* Data type for storing 2D greyscale image */
typedef struct imagestruct
{
  int width, height;
  int **imdata;
} *Image;

static void error(char *errmsg)
{ /* print error message an abort program */
  fprintf (stderr, errmsg);
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
  fscanf (f, "%d %d\n", &w, &h);
  /* read maximum greyvalue (dummy) */
  fscanf (f, "%d\n", &maxval);
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
    fread(scanline, 1, w, f);
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

static void contrastStretch (int low, int high, Image image)
{
  int row, col, min, max;
  int width=image->width, height=image->height, **im=image->imdata;
  float scale;
    
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
  scale = (float)(high-low) / (max-min);

  /* Stretch image */
  for (row = 0; row < height; row++)
  {
    for (col = 0; col < width; col++)
    {
      im[row][col] = scale*(im[row][col] - min);
    }
  }
}


int main(int argc, char **argv)
{
  Image image;
  
  /* Initialize MPI */
  MPI_Init (&argc, &argv);

  if (argc != 3)
  {
    printf ("Usage: %s input.pgm output.pgm\n", argv[0]);
    exit (EXIT_FAILURE);
  }

  image = readImage(argv[1]);
  contrastStretch(0, 255, image);
  writeImage(image, argv[2]);

  freeImage(image);
  MPI_Finalize ();
  return EXIT_SUCCESS;
}
