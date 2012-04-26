/*
 * repmat.c
 *
 * Code generation for function 'repmat'
 *
 * C source code generated on: Tue Apr 24 12:18:22 2012
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "convertMovieToBinaryArray.h"
#include "repmat.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */
static emlrtBCInfo d_emlrtBCI = { 1, 275558400, 64, 21, "", "repmat", "/Applications/MATLAB_R2011b.app/toolbox/eml/lib/matlab/elmat/repmat.m", 0 };
static emlrtBCInfo e_emlrtBCI = { 1, 307200, 64, 29, "", "repmat", "/Applications/MATLAB_R2011b.app/toolbox/eml/lib/matlab/elmat/repmat.m", 0 };

/* Function Declarations */

/* Function Definitions */

void repmat(const uint8_T a[307200], uint8_T b[275558400])
{
    int32_T ib;
    int32_T jtilecol;
    int32_T iacol;
    int32_T jcol;
    int32_T k;
    ib = 1;
    for (jtilecol = 0; jtilecol < 897; jtilecol++) {
        iacol = 1;
        for (jcol = 0; jcol < 640; jcol++) {
            for (k = 0; k < 480; k++) {
                b[emlrtBoundsCheck(ib, &d_emlrtBCI) - 1] = a[emlrtBoundsCheck(iacol, &e_emlrtBCI) - 1];
                iacol++;
                ib++;
            }
        }
    }
}
/* End of code generation (repmat.c) */
