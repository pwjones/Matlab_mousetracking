/*
 * mean.c
 *
 * Code generation for function 'mean'
 *
 * C source code generated on: Tue Apr 24 12:18:22 2012
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "convertMovieToBinaryArray.h"
#include "mean.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */
static emlrtRSInfo g_emlrtRSI = { 36, "mean", "/Applications/MATLAB_R2011b.app/toolbox/eml/lib/matlab/datafun/mean.m" };
static emlrtBCInfo emlrtBCI = { 1, 307200, 75, 13, "", "sum", "/Applications/MATLAB_R2011b.app/toolbox/eml/lib/matlab/datafun/sum.m", 0 };
static emlrtBCInfo b_emlrtBCI = { 1, 275558400, 69, 22, "", "sum", "/Applications/MATLAB_R2011b.app/toolbox/eml/lib/matlab/datafun/sum.m", 0 };
static emlrtBCInfo c_emlrtBCI = { 1, 275558400, 72, 30, "", "sum", "/Applications/MATLAB_R2011b.app/toolbox/eml/lib/matlab/datafun/sum.m", 0 };

/* Function Declarations */

/* Function Definitions */

void mean(const uint8_T x[275558400], real_T y[307200])
{
    int32_T iy;
    int32_T ixstart;
    int32_T j;
    int32_T ix;
    real_T s;
    int32_T k;
    EMLRTPUSHRTSTACK(&g_emlrtRSI);
    iy = 0;
    ixstart = 0;
    for (j = 0; j < 307200; j++) {
        ixstart++;
        ix = ixstart;
        s = (real_T)x[emlrtBoundsCheck(ixstart, &b_emlrtBCI) - 1];
        for (k = 0; k < 896; k++) {
            ix += 307200;
            s += (real_T)x[emlrtBoundsCheck(ix, &c_emlrtBCI) - 1];
        }
        iy++;
        y[emlrtBoundsCheck(iy, &emlrtBCI) - 1] = s;
    }
    for (iy = 0; iy < 307200; iy++) {
        y[iy] /= 897.0;
    }
    EMLRTPOPRTSTACK(&g_emlrtRSI);
}
/* End of code generation (mean.c) */
