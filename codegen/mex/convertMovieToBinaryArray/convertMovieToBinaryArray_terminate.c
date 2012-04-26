/*
 * convertMovieToBinaryArray_terminate.c
 *
 * Code generation for function 'convertMovieToBinaryArray_terminate'
 *
 * C source code generated on: Tue Apr 24 12:18:22 2012
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "convertMovieToBinaryArray.h"
#include "convertMovieToBinaryArray_terminate.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */

/* Function Definitions */

void convertMovieToBinaryArray_atexit(void)
{
    emlrtExitTimeCleanup(&emlrtContextGlobal);
}

void convertMovieToBinaryArray_terminate(void)
{
    emlrtLeaveRtStack(&emlrtContextGlobal);
}
/* End of code generation (convertMovieToBinaryArray_terminate.c) */
