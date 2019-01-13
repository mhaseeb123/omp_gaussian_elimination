/*
 * Author: Muhammad Haseeb
 * Class: CS6260 | Spring-2018
 * File: GaussianElimination.cpp
 * Assignment: Gaussian Elimination
 * Date: 03/15/2018
 *
 */

/* Include file */
#include "GaussianElimination.h"

#ifdef DEBUG
 #ifdef WINDOWS
  #include <Windows.h>
 #endif /* WINDOWS */
#endif /* DEBUG */

using namespace std;
/*
 * FUNCTION: main
 * DESCRIPTION:
 *
 * @params: argc, argv
 * @returns: status
 */
STATUS main(int argc, char** argv)
{
    STATUS status = SUCCESS;

    /* The size of A matrix */
    int n;
    /* Used for indexing into arrays */
    int i = 0;
    int j = 0;
    int k = 0;
    int row = 0;
    int col = 0;

    /* Factor calculation in the Echelon reduction */
    double factor = 0.0;

    /* Times */
    double start_time = 0.0;
    double gauss_end = 0.0;
    double subs_start = 0.0;
    double subs_end = 0.0;

    /* Check if number of tosses passed as parameter */
    if (argc < 3)
    {

        cout << "ERROR: Number of threads and/or Matrix dimension not specified\n";
        cout << "USAGE: ./GaussianElimination.exe <NUM_THREADS> <DIMENSION>\n";

        status = ERR_INVLD_ARGS;
        exit(0);
    }

    UINT thread_count = atoi(argv[1]);
    n = atoi(argv[2]);

    if (n<=0)
    {
        n = DEFAULT_N;
    }

    if (thread_count <=0)
    {
        thread_count = DEFAULT_THREADS;
    }

    /* Allocate memory for the n x n matrix system */
    double *A = new double[n * n];
    double *b = new double[n];
    double *x = new double[n];

    /* TODO: Status and sanity checks */

    /* Initialize the system */
    Init(A, b, x, n);

    cout << "Init Done:" << endl;
    cout << "Threads: " << thread_count << endl;
    cout << "Dimension: "<< n << endl;

#ifdef DEBUG
    PrintMatrix(A, n);
    PrintVector(x, n, 'x');
    PrintVector(b, n, 'b');
#endif /* DEBUG */

/*************************** GAUSSIAN *********************************/
    /* Mark the start of Gaussian Elimination here */
    start_time = omp_get_wtime();

    /* Gaussian Elimination code here */
    for (i = 0; i < n - 1; i++)
    {
        /* HM: This loop can be parallelized.
         *     Not dependent on anything from
         *     any other iteration
         */
#pragma omp parallel for num_threads(thread_count) schedule(runtime)  \
    default(none) private(j, k, factor) shared(i,A,b,n,x,thread_count)
        for (j = i + 1; j < n; j++)
        {
            factor = A[j * n + i] / A[i * n + i];

            for (k = 0; k < n; k++)
            {
                if (k <= i)
                {
                    A[j * n + k] = 0;
                }
                else
                {
                    A[j * n + k] -= factor * A[i * n + k];
                }
            }

            b[j] -= factor * b[i];
        }
    }

#ifdef DEBUG
    PrintMatrix(A, n);
    PrintVector(x, n, 'x');
    PrintVector(b, n, 'b');
#endif /* DEBUG */

    gauss_end = omp_get_wtime();
    cout << "Gaussian Elimination Done: " << endl;
    cout << "Elapsed Time: " << (float)(gauss_end - start_time) << 's' << endl << endl;

/********************* BACK SUBSTITUTION ********************************/
    subs_start = omp_get_wtime();

    /* Row-oriented back substitution code here */
    for (row = n - 1; row >= 0; row--)
    {
        x[row] = b[row];

#pragma omp parallel for num_threads(thread_count) schedule(runtime)  \
    default(none) private(col) shared(A,b,n,x,row)
        for (col = row + 1; col < n; col++)
        {
#pragma omp atomic
            x[row] -= A[row * n + col] * x[col];
        }

        x[row] /= A[row * n + row];
    }

    subs_end = omp_get_wtime();
    cout << "Back Substitution Done: " << endl;
    cout << "Elapsed Time: " << (float)(subs_end - subs_start) << 's' << endl << endl;

    #ifdef DEBUG
    PrintMatrix(A, n);
    PrintVector(x, n, 'x');
    PrintVector(b, n, 'b');
#endif /* DEBUG */

/********************** CLEANUP ***********************/
    /* Deallocate memory */
    delete[] A;
    delete[] x;
    delete[] b;

    /* Set pointers to NULL */
    A = NULL;
    x = NULL;
    b = NULL;

    cout << "All done" << endl;
    cout << "Total Elapsed Time: " << (float)(subs_end - start_time) << "s\n\n";
    return status;
}

/*
 * FUNCTION: Init
 * DESCRIPTION: Assume that A[] is n*n matrix, b[] & x[] are vectors of n elements.
 *
 * @params: argc, argv
 * @returns: status
 */
void Init(double *A, double *b, double *x, int n)
{
    int i;
    int j;

    for (i = 0; i < n; i++)
    {
        x[i] = 1.0;
    }

    srand(1);

    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
            if (i != j)
                A[i * n + j] = rand() / ((double) RAND_MAX);
            else
                A[i * n + i] = n / 10.0;
    }

    for (i = 0; i < n; i++)
    {
        b[i] = 0;
        for (j = 0; j < n; j++)
            b[i] += A[i * n + j] * x[j];
    }

    memset((void *)x, 0x0, n * sizeof(double));
}

/* Funtion to Print a matrix */
void PrintMatrix(double *A, int n)
{
    int i = 0;
    int j = 0;

    cout << endl << "Printing Matrix A:"<<endl;

    /* Print the matrix here */
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            cout << A[i * n + j] << '\t';
        }

        cout << endl;
    }
}

/* Funtion to Print a vector */
void PrintVector(double *A, int n, char nm)
{
    int i = 0;

    cout << endl << "Printing Vector: "<< nm <<endl;

    /* Print the matrix here */
    for (i = 0; i < n; i++)
    {
        cout << A[i] << endl;
    }
}

/*
 * FUNCTION: Sample Init
 * DESCRIPTION: Assume that A[] is 3*3 matrix, b[] & x[] are vectors of n elements.
 *
 * @params: A, b, x, n
 * @returns: status
 */
void Sample_Init(double *A, double *b, double *x)
{
    A[0] = 2;
    A[1] = -3;
    A[2] = 0;
    A[3] = 4;
    A[4] = -5;
    A[5] = 1;
    A[6] = 2;
    A[7] = -1;
    A[8] = -3;

    x[0] = 1;
    x[1] = 1;
    x[2] = 1;

    b[0] = 3;
    b[1] = 7;
    b[2] = 5;
}
