/*
 * convertMovieToBinaryArray_mexutil.c
 *
 * Code generation for function 'convertMovieToBinaryArray_mexutil'
 *
 * C source code generated on: Tue Apr 24 12:18:22 2012
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "convertMovieToBinaryArray.h"
#include "convertMovieToBinaryArray_mexutil.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */

/* Function Definitions */

const mxArray *b_emlrt_marshallOut(const uint8_T u[275558400])
{
    const mxArray *y;
    static const int32_T iv1[3] = { 480, 640, 897 };
    const mxArray *m1;
    uint8_T (*pData)[];
    int32_T i;
    y = NULL;
    m1 = mxCreateNumericArray(3, (int32_T *)&iv1, mxUINT8_CLASS, mxREAL);
    pData = (uint8_T (*)[])mxGetData(m1);
    for (i = 0; i < 275558400; i++) {
        (*pData)[i] = u[i];
    }
    emlrtAssign(&y, m1);
    return y;
}
/* End of code generation (convertMovieToBinaryArray_mexutil.c) */
