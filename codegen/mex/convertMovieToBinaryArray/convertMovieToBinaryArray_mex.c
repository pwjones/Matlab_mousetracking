/*
 * convertMovieToBinaryArray_mex.c
 *
 * Code generation for function 'convertMovieToBinaryArray'
 *
 * C source code generated on: Tue Apr 24 12:18:22 2012
 *
 */

/* Include files */
#include "mex.h"
#include "convertMovieToBinaryArray_api.h"
#include "convertMovieToBinaryArray_initialize.h"
#include "convertMovieToBinaryArray_terminate.h"

/* Type Definitions */

/* Function Declarations */
static void convertMovieToBinaryArray_mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);

/* Variable Definitions */
emlrtContext emlrtContextGlobal = { true, false, EMLRT_VERSION_INFO, NULL, "convertMovieToBinaryArray", NULL, false, NULL, false, 1, false };

/* Function Definitions */
static void convertMovieToBinaryArray_mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  /* Temporary copy for mex outputs. */
  mxArray *outputs[1];
  int n = 0;
  int nOutputs = (nlhs < 1 ? 1 : nlhs);
  c_convertMovieToBinaryArrayStac* d_convertMovieToBinaryArrayStac = (c_convertMovieToBinaryArrayStac*)mxCalloc(1,sizeof(c_convertMovieToBinaryArrayStac));
  /* Check for proper number of arguments. */
  if(nrhs != 1) {
    mexErrMsgIdAndTxt("emlcoder:emlmex:WrongNumberOfInputs","1 input required for entry-point 'convertMovieToBinaryArray'.");
  } else if(nlhs > 1) {
    mexErrMsgIdAndTxt("emlcoder:emlmex:TooManyOutputArguments","Too many output arguments for entry-point 'convertMovieToBinaryArray'.");
  }
  /* Module initialization. */
  convertMovieToBinaryArray_initialize(&emlrtContextGlobal);
  /* Call the function. */
  convertMovieToBinaryArray_api(d_convertMovieToBinaryArrayStac, prhs,(const mxArray**)outputs);
  /* Copy over outputs to the caller. */
  for (n = 0; n < nOutputs; ++n) {
    plhs[n] = emlrtReturnArrayR2009a(outputs[n]);
  }
  /* Module finalization. */
  convertMovieToBinaryArray_terminate();
  mxFree(d_convertMovieToBinaryArrayStac);
}
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  /* Initialize the memory manager. */
  mexAtExit(convertMovieToBinaryArray_atexit);
  emlrtClearAllocCount(&emlrtContextGlobal, 0, 0, NULL);
  /* Dispatch the entry-point. */
  convertMovieToBinaryArray_mexFunction(nlhs, plhs, nrhs, prhs);
}
/* End of code generation (convertMovieToBinaryArray_mex.c) */
