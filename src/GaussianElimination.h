/*
 * Author: Muhammad Haseeb
 * Class: CS6260 | Spring-2018
 * File: GaussianElimination.h
 * Assignment: Gaussian Elimination
 * Date: 03/15/2018
 *
 */

#ifndef GAUSSIANELIMINATION_H_
#define GAUSSIANELIMINATION_H_

/* Include STD files */
#include <iostream>
#include <cstring>
#include <cstdlib>
#include <ctime>
#include <iomanip>

using namespace std;

/* Include OpenMP */
#ifdef _OPENMP
 #include <omp.h>
#endif /* _OPENMP */

/* Enable if debug info required - Works with Windows only */
#undef DEBUG
#define WINDOWS

/* Error codes */
#define SUCCESS                0
#define ERR_BASE               -3000
#define ERR_PROC_ID            (ERR_BASE - 1)
#define ERR_INVLD_ARGS         (ERR_BASE - 2)

/* Defines */
#ifdef _OPENMP
 #define DEFAULT_THREADS           8
#endif /* _OPENMP */

#define DEFAULT_N                  1000

/* Macro to avoid compiler warnings if a parameter is not being used */
#define UNUSED_PARAM(arg)   (void) arg

/* Typedefs */
typedef unsigned int   UINT;
typedef unsigned char  UCHAR;
typedef unsigned short USHORT;
typedef int            STATUS;


/* Function Prototypes */
STATUS main(int argc, char** argv);
void Init(double *A, double *b, double *x, int n);
void PrintMatrix(double *A, int n);
void PrintVector(double *A, int n, char nm);
void Sample_Init(double *A, double *b, double *x);

#endif /* GAUSSIANELIMINATION_H_ */
