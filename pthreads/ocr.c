/* file: ocr.c
 *(C) Arnold Meijster and Rob de Bruin
 *
 * A simple OCR demo program.
 */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define FALSE 0
#define TRUE  1

#define FOREGROUND 255
#define BACKGROUND 0

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

#define CORRTHRESHOLD 0.95

/* Data type for storing 2D greyscale image */
typedef struct imagestruct
{
    int width, height;
    int **imdata;
} *Image;

int charwidth, charheight;  /* width an height of templates/masks */

#define NSYMS 75
Image mask[NSYMS];          /* templates: one for each symbol     */
char  symbols[NSYMS] = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
"abcdefghijklmnopqrstuvwxyz"
"0123456789-+*/,.?!:;'()";

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

static void threshold(int th, int less, int greater, Image image)
{ /* routine for converting a greyscale image into a binary image
   * (in place) by means of thresholding. The first parameter is
   * the threshold value, the  second parameter is the value in the
   * output image for pixels that were less than the threshold,
   * while the third parameter is the value in the output for pixels
   * that were at least the threshold value.
   */
    int row, col;
    int width=image->width, height=image->height, **im=image->imdata;
    for (row = 0; row < height; row++)
    {
        for (col = 0; col < width; col++)
        {
            im[row][col] = (im[row][col] < th ? less : greater);
        }
    }
}

static void distanceBlur(int background, Image image)
{ /* blur the image by assigning to each pixel a gray value
   * which is inverse proportional to the distance to
   * its nearest foreground pixel
   */
    int width=image->width, height=image->height, **im=image->imdata;
    int r, c, dt;
    /* forward pass */
    im[0][0] = (im[0][0] == background ? width + height : 0);
    for (c = 1 ; c < width; c++)
    {
        im[0][c] = (im[0][c] == background ? 1 + im[0][c-1] : 0);
    }
    /* other scanlines */
    for (r = 1; r < height; r++)
    {
        im[r][0] = (im[r][0] == background ? 1 + im[r-1][0] : 0);
        for (c = 1; c < width; c++)
        {
            im[r][c] = (im[r][c] == background ?
                    1 + MIN(im[r][c-1], im[r-1][c]) : 0);
        }
    }
    /* backward pass */
    dt = im[height-1][width-1];
    for (c = width-2; c >= 0; c--)
    {
        im[height-1][c] = MIN(im[height-1][c], im[height-1][c+1]+1);
        dt = MAX(dt, im[height-1][c]);
    }
    /* other scanlines */
    for (r = height-2; r >= 0; r--)
    {
        im[r][width-1] = MIN(im[r][width-1], im[r+1][width-1]+1);
        dt = MAX(dt, im[r][width-1]);
        for (c = width-2; c >= 0; c--)
        {
            im[r][c] = MIN(im[r][c], 1 + MIN(im[r+1][c], im[r][c+1]));
            dt = MAX(dt, im[r][c]);
        }
    }
    /* At this point the grayvalue of each pixel is
     * the distance to its nearest foreground pixel in the original
     * image. Also, dt is the maximum distance obtained.
     * We now 'invert' this image, by making the distance of the
     * foreground pixels maximal.
     */
    for (r = 0; r < height; r++)
    {
        for (c = 0; c < width; c++)
        {
            im[r][c] = dt - im[r][c];
        }
    }
}

static double PearsonCorrelation(Image image, Image mask)
{ /* returns the Pearson correlation coefficient of image and mask. */
    int width=image->width, height=image->height;
    int r, c, **x=image->imdata, **y=mask->imdata;
    double mx=0, my=0, sx=0, sy=0, sxy=0;
    /* Calculate the means of x and y */
    for (r = 0; r < height; r++)
    {
        for (c = 0; c < width; c++)
        {
            mx += x[r][c];
            my += y[r][c];
        }
    }
    mx /= width*height;
    my /= width*height;
    /* Calculate correlation */
    for (r = 0; r < height; r++)
    {
        for (c = 0; c < width; c++)
        {
            sxy += (x[r][c] - mx)*(y[r][c] - my);
            sx  += (x[r][c] - mx)*(x[r][c] - mx);
            sy  += (y[r][c] - my)*(y[r][c] - my);
        }
    }
    return (sxy / sqrt(sx*sy));
}

static int PearsonCorrelator(Image character)
{
    /* returns the index of the best matching symbol in the
     * alphabet array. A negative value means that the correlator
     * is 'uncertain' about the match, but the absolute value is still
     * considered to be the best match.
     */
    double best = PearsonCorrelation(character, mask[0]);
    int sym, bestsym = 0;
    for (sym=1; sym<NSYMS; sym++)
    {
        double corr = PearsonCorrelation(character, mask[sym]);
        if (corr > best)
        {
            best = corr;
            bestsym = sym;
        }
    }
    return (best >= CORRTHRESHOLD ? bestsym : -bestsym);
}

static int isEmptyRow(int row, int col0, int col1,
        int background, Image image)
{ /* returns TRUE if and only if all pixels im[row][c]
   * (with col0<=c<col1) are background pixels
   */
    int col, **im=image->imdata;
    for (col = col0; col < col1; col++)
    {
        if (im[row][col] != background)
        {
            return FALSE;
        }
    }
    return TRUE;
}

static int isEmptyColumn(int row0, int row1, int col,
        int background, Image image)
{ /* returns TRUE if and only if all pixels im[r[col]
   * (with row0<=r<row1) are background pixels
   */
    int row, **im=image->imdata;
    for (row = row0; row < row1; row++)
    {
        if (im[row][col] != background)
        {
            return FALSE;
        }
    }
    return TRUE;
}

static void makeCharImage(int row0, int col0, int row1, int col1,
        int background, Image image, Image mask)
{ /* construct a (centered) mask image from a segmented character */
    int r, r0, r1, c, c0, c1, h, w;
    int **im = image->imdata, **msk = mask->imdata;
    /* When a character has been segmented, its left and right margin
     * (given by col1 and col1) are tight, however the top and bottom
     * margin need not be. This routine tightens the top and
     * bottom marging as well.
     */
    while (isEmptyRow(row0, col0, col1, background, image))
    {
        row0++;
    }
    while (isEmptyRow(row1-1, col0, col1, background, image))
    {
        row1--;
    }
    /* Copy image data (centered) into a mask image */
    h = row1 - row0;
    w = col1 - col0;
    r0 = MAX(0, (charheight-h)/2);
    r1 = MIN(charheight, r0 + h);
    c0 = MAX(0, (charwidth - w)/2);
    c1 = MIN(charwidth, c0 + w);
    for (r=0; r<charheight; r++)
    {
        for (c=0; c<charwidth; c++)
        {
            msk[r][c] = ((r0 <= r) && (r < r1) && (c0 <= c) && (c < c1) ?
                    im[r-r0+row0][c-c0+col0] : background);
        }
    }
    /* blur mask */
    distanceBlur(background, mask);
}

static int findCharacter(int background, Image image,
        int row0, int row1, int *col0, int *col1)
{ /* find the bounding box of the left most character
   * with column>=col0 in the linestrip given by row0 and row 1.
   * Note that col0 is an input/output parameter, while
   * col1 is a strict output parameter.
   * This routine returns TRUE if a character was found,
   * and FALSE otherwise (at end of line).
   */

    /* find fist column which has at least one foreground (ink) pixel */
    int width=image->width;
    while ((*col0 < width) &&
            (isEmptyColumn(row0, row1, *col0, background, image)))
    {
        (*col0)++;
    }
    /* find column which is entirely empty (paper) */
    *col1 = *col0;
    while ((*col1 < width) &&
            (!isEmptyColumn(row0, row1, *col1, background, image)))
    {
        (*col1)++;
    }
    return (*col0 < *col1 ? TRUE : FALSE);
}

static void characterSegmentation(int background,
        int row0, int row1, Image image)
{ /* Segment a linestrip into characters and perform
   * character matching using a Pearson correlator.
   * This routine also prints the output characters.
   */
    int match, col1, prev=0, col0=0;
    Image token = makeImage(charwidth, charheight);
    while (findCharacter(background, image, row0, row1, &col0, &col1))
    {
        /* Was there a space ? */
        if (prev > 0)
        {
            int i;
            for (i = 0; i < (int)((col0-prev)/(1.1*charwidth)); i++)
            {
                printf (" ");
            }
        }
        /* match character */
        makeCharImage(row0, col0, row1, col1, background, image, token);
        match = PearsonCorrelator(token);
        if (match >= 0)
        {
            printf ("%c", symbols[match]);
        } else
        {
            printf ("#%c#", symbols[-match]);
        }
#if 1
        {
            /* Mark bounding box of character.
             * This can be turned on for debugging pruposes.
             * In combination with writePGM it is possible
             * to write the 'segmnented' image to a file.
             */
            int **im=image->imdata;
            int r, c;
            for (r = row0; r < row1; r++)
            {
                for (c = col0; c < col1; c++)
                {
                    if (im[r][c] == background)
                    {
                        im[r][c] = 128;
                    }
                }
            }
        }
#endif
        col0 = prev = col1;
    }
    freeImage(token);
    return;
}

static int findLineStrip (int background, int *row0, int *row1,
        Image image)
{ /* find the first line strip that is encountered
   * when parsing scanlines from top to bottom starting
   * from row0. Note that row0 is an input/output parameter, while
   * row1 is a strict output parameter.
   * This routine returns TRUE if a line strip was found,
   * and FALSE otherwise (at end of page).
   */
    int width=image->width, height=image->height;
    /* find scanline which has at least one foreground (ink) pixel */
    while ((*row0 < height) &&
            (isEmptyRow(*row0, 0, width, background, image)))
    {
        (*row0)++;
    }
    /* find scanline which is entirely empty (paper) */
    *row1 = *row0;
    while ((*row1 < height) &&
            (!isEmptyRow(*row1, 0, width, background, image)))
    {
        (*row1)++;
    }
    return (*row0 < *row1 ? TRUE : FALSE);
}

static void lineSegmentation(int background, Image page)
{ /* Segments a page into line strips. For each line strip
   * the character recognition pipeline is started.
   */
    int row1, row0=0, prev=0;
    while (findLineStrip(background, &row0, &row1, page))
    {
        /* Was there an empty line? */
        if (prev > 0)
        {
            int i;
            for (i = 0; i < (int)((row0 - prev)/(1.2*charheight)); i++)
            {
                printf ("\n");
            }
        }
        /* separate characters in line strip */
        characterSegmentation(background, row0, row1, page);
        row0 = prev = row1;
        printf ("\n");
    }
#if 0
    /* You can enable this code fragment for debugging purposes */
    writePGM(page, "segmentation.pgm");
#endif
}

static void constructAlphabetMasks()
{ /* Construct a full set of alphabet symbols from the file
   * 'alphabet.pgm'. The symbols are represented in the
   * template images mask[].
   */
    int row0, row1, col0, col1, sym;
    Image alphabet = readPGM("alphabet.pgm");
    threshold(100, FOREGROUND, BACKGROUND, alphabet);

    /* first pass: only used for computing bounding box of largest
     * character in the alphabet. */
    charwidth = charheight = 0;
    row0 = sym = 0;
    while (sym < NSYMS)
    {
        if (!findLineStrip (BACKGROUND, &row0, &row1, alphabet))
        {
            error("Error: construction of alphabet symbols failed.\n");
        }
        charheight = MAX(charheight, row1 - row0);
        col0 = 0;
        while ((sym < NSYMS) &&
                (findCharacter(BACKGROUND, alphabet, row0, row1, &col0, &col1)))
        {
            /* we found a token with bounding box: (row0,col0)-(row1,col1) */
            charwidth = MAX(charwidth, col1-col0);
            sym++;
            col0 = col1;
        }
        row0 = row1;
    }

    /* allocate space for masks */
    for (sym=0; sym<NSYMS; sym++)
    {
        mask[sym] = makeImage(charwidth, charheight);
    }
    /* second pass */
    row0 = sym = 0;
    while (sym < NSYMS)
    {
        findLineStrip (BACKGROUND, &row0, &row1, alphabet);
        col0 = 0;
        while ((sym < NSYMS) &&
                (findCharacter(BACKGROUND, alphabet, row0, row1, &col0, &col1)))
        {
            makeCharImage(row0, col0, row1, col1,
                    BACKGROUND, alphabet, mask[sym]);
#if 1
            {
                /* You can enable this code fragment for debugging purposes */
                char filename[16];
                sprintf(filename, "template%02d.pgm", sym);
                writePGM(mask[sym], filename);
            }
#endif
            sym++;
            col0 = col1;
        }
        row0 = row1;
    }
    freeImage(alphabet);
}

int main(int argc, char **argv)
{
    int i;
    Image image;

    if (argc < 2)
    {
        fprintf(stderr, "Usage: %s page1.pgm page2.pgm...\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    /* process alphabet image */
    constructAlphabetMasks();

    /* process pages */
    for (i = 1; i < argc; i++)
    {
        printf ("\n--%s----------------------------------\n", argv[i]);
        image = readPGM(argv[i]);
        threshold(100, FOREGROUND, BACKGROUND, image);
        lineSegmentation(BACKGROUND, image);
        freeImage(image);
    }

    /* clean up memory used for alphabet masks */
    for (i = 0; i < NSYMS; i++)
    {
        freeImage(mask[i]);
    }
    return EXIT_SUCCESS;
}
