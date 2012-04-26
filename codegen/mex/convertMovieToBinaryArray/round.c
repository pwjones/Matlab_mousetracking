/*
 * round.c
 *
 * Code generation for function 'round'
 *
 * C source code generated on: Tue Apr 24 12:18:22 2012
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "convertMovieToBinaryArray.h"
#include "round.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */

/* Function Definitions */

void b_round(real_T x[307200])
{
    int32_T k;
    real_T b_x;
    for (k = 0; k < 307200; k++) {
        if (muDoubleScalarAbs(x[k]) > 4.503599627370496E+15) {
            b_x = x[k];
        } else if (x[k] >= 0.5) {
            b_x = muDoubleScalarFloor(x[k] + 0.5);
        } else if (x[k] > -0.5) {
            b_x = x[k] * 0.0;
        } else {
            b_x = muDoubleScalarCeil(x[k] - 0.5);
        }
        x[k] = b_x;
    }
}
/* End of code generation (round.c) */
