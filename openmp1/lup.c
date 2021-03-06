/* File: lup.c
 * Author: Arnold Meijster (a.meijster@rug.nl)
 * Version: 1.0 (01 March 2008)
 */
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <omp.h>
#include <string.h>

#define real float
#define ABS(a) ((a)<0 ? (-(a)) : (a))

/* forward references */
static void *safeMalloc(size_t n);
static real **allocMatrix(size_t height, size_t width);
static void freeMatrix(real **mat);
static real *allocVector(size_t length);
static void freeVector(real *vec);
static void showMatrix (size_t n, real **A);
static void showVector (size_t n, real *vec);
static void decomposeLUP(size_t n, real **A, size_t *P);
static void LUPsolve(size_t n, real **LU, size_t *P, real *x, real *b);
static void solve(size_t n, real **A, real *x, real *b);
static void invert(size_t n, real**A);
static real** multiply(size_t n, real **A, real **B);
/********************/

void *safeMalloc(size_t n)
{
    void *ptr;
    ptr = malloc(n);
    if (ptr == NULL)
    {
        fprintf (stderr, "Error: malloc(%lu) failed\n", n);
        exit(-1);
    }
    return ptr;
}

real **allocMatrix(size_t height, size_t width)
{
    real **matrix;
    size_t row;

    matrix = safeMalloc(height * sizeof(real *));
    matrix[0] = safeMalloc(width*height*sizeof(real));
    for (row=1; row<height; ++row)
        matrix[row] = matrix[row-1] + width;
    return matrix;
}

void freeMatrix(real **mat)
{
    free(mat[0]);
    free(mat);
}

real *allocVector(size_t length)
{
    return safeMalloc(length*sizeof(real));
}


void freeVector(real *vec)
{
    free(vec);
}


void showMatrix (size_t n, real **A)
{
    size_t i, j;
    for (i=0; i<n; ++i)
    {
        for (j=0; j<n; ++j)
        {
            printf ("%f ", A[i][j]);
        }
        printf ("\n");
    }
}

void showVector (size_t n, real *vec)
{
    size_t i;
    for (i=0; i<n; ++i)
    {
        printf ("%f ", vec[i]);
    }
    printf ("\n");
}

void decomposeLUP(size_t n, real **A, size_t *P)
{
    /* computes L, U, P such that A=L*U*P */

    int h, i, j, k, row;
    real pivot, absval, tmp;

    #pragma omp parallel for if (n>100)
    for (i=0; i<n; ++i)
    {
        P[i] = i;
    }

    //#pragma omp parallel for private(h,i,j,k,row,pivot,absval,tmp) schedule(static,1)
    for (k=0; k<n-1; ++k)
    {
        row = -1;
        pivot = 0;

        for (i=k; i<n; ++i)
        {
            //absval = (A[i][k] >= 0 ? A[i][k] : -A[i][k]);
            absval = ABS(A[i][k]);
            if (absval>pivot)
            {
                pivot = absval;
                row = i;
            }
        }

        if (row == -1)
        {
            printf ("Singular matrix\n");
            exit(-1);
        }

        /* swap(P[k],P[row]) */
        h = P[k];
        P[k] = P[row];
        P[row] = h;

        /* swap rows */
        #pragma omp parallel for private(tmp,i) if (n>100)
        for (i=0; i<n; ++i)
        {
            tmp = A[k][i];
            A[k][i] = A[row][i];
            A[row][i] = tmp;
        }

        #pragma omp parallel for private(i,j)
        for (i=k+1; i<n; ++i)
        {
            A[i][k] /= A[k][k];
            for (j=k+1; j<n; ++j)
            {
                A[i][j] -= A[i][k]*A[k][j];
            }
        }
    }
}

void LUPsolve(size_t n, real **LU, size_t *P, real *x, real *b)
{
    real *y, s;
    size_t i, j;

    /* Solve Ly=Pb using forward substitution */
    y = x;  /* warning, y is an alias for x! It is safe, though. */
    for (i=0; i<n; ++i)
    {
        s = 0;
        y[i] = b[P[i]];
        #pragma omp parallel for reduction(+:s)
        for (j=0; j<i; ++j)
        {
            s += LU[i][j]*y[j];
        }
        y[i] -= s;
    }

    /* Solve Ux=y using backward substitution */
    i=n;
    while (i>0)
    {
        i--;
        x[i] = y[i];
        s = 0;
        #pragma omp parallel for reduction(+:s)
        for (j=i+1; j<n; ++j)
        {
            s += LU[i][j]*x[j];
        }
        x[i] -= s;
        x[i] /= LU[i][i];
    }
}

void solve(size_t n, real **A, real *x, real *b)
{
    size_t *P;

    /* Construct LUP decomposition */
    P = safeMalloc(n*sizeof(size_t));
    decomposeLUP(n, A, P);
    /* Solve by forward and backward substitution */
    LUPsolve(n, A, P, x, b);

    free(P);
}

void invert(size_t n, real **A)
{
    int i, k;
    real **Ainv, *b;
    size_t *P;

    Ainv = allocMatrix(n, n);
    P = safeMalloc(n*sizeof(size_t));
    /* Start by constructing a LUP decomposition */
    decomposeLUP(n, A, P);
    /* Invert matrix by solving for each column of the identity matrix */

    #pragma omp parallel private(b)
    {
    b = allocVector(n);
    #pragma omp for
    for (k=0; k<n; ++k)
    {
        for (i=0; i<n; ++i)
            b[i] = 0;
        b[k] = 1;
        LUPsolve(n, A, P, Ainv[k], b);
    }
    freeVector(b);
    }

    #pragma omp parallel for private(i,k) schedule(static,1)
    for (i=0; i<n; ++i)
    {
        for (k=0; k<n; ++k)
            A[i][k] = Ainv[k][i];
    }

    freeMatrix(Ainv);
    //freeVector(b);
    free(P);
}

real** multiply(size_t n, real **A, real **B)
{
    int i, j, k;
    real **C, s;

    C = allocMatrix(n, n);

    #pragma omp parallel for private(i,j,k,s) schedule(static, 1)
    for (i=0; i<n; ++i)
    {
        #pragma omp parallel for private(j,k,s) schedule(static, 1)
        for (j=0; j<n; ++j)
        {
            s = 0;
            #pragma omp parallel for reduction(+:s)
            for (k=0; k<n; ++k)
            {
                s += A[i][k] * B[k][j];
            }
            C[i][j] = s;
        }
    }
    return C;
}

int main(int argc, char **argv)
{
    int n, i, j;
    real **A, **B, *x, *b;

    if (argc == 2) {
        n = atoi(argv[1]);
        A = allocMatrix(n, n);
        x = allocVector(n);
        b = allocVector(n);

        #pragma omp parallel for private(i,j) schedule(static, 1)
        for(i=0; i<n; ++i) {
            b[i] = 1;
            #pragma omp parallel for private(j) schedule(static, 1)
            for (j=0; j<n; ++j) {
                if (i == j)
                    A[i][j] = -2;
                else
                    A[i][j] = (ABS(i-j) == 1);
            }
        }
    } else {
        n = 3;
        A = allocMatrix(n, n);
        x = allocVector(n);
        b = allocVector(n);

        A[0][0]= 0; A[0][1]= 2; A[0][2]=  3;
        A[1][0]= 3; A[1][1]= 0; A[1][2]=  1;
        A[2][0]= 6; A[2][1]= 2; A[2][2]=  8;

        b[0] = 5;
        b[1] = 4;
        b[2] = 16;
    }
#ifdef INVERT
    B = allocMatrix(n, n);
    for (i=0; i<n; ++i)
        for (j=0; j<n; ++j)
            B[i][j] = A[i][j];
#endif

    showMatrix(3, A);
    printf ("\n");

#ifndef INVERT
    solve(n, A, x, b);

    showVector(3, x);
#else
    invert(n, A);
    showMatrix(3, B);
    showMatrix(3, A);
    printf("\n");

    showMatrix(10, multiply(n, A, B));
#endif

    freeMatrix(A);
    freeVector(x);
    freeVector(b);

    return EXIT_SUCCESS;
}
