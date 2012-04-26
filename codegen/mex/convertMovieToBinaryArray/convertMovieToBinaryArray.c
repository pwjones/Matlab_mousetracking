/*
 * convertMovieToBinaryArray.c
 *
 * Code generation for function 'convertMovieToBinaryArray'
 *
 * C source code generated on: Tue Apr 24 12:18:22 2012
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "convertMovieToBinaryArray.h"
#include "repmat.h"
#include "round.h"
#include "mean.h"
#include "convertMovieToBinaryArray_mexutil.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */
static emlrtRSInfo emlrtRSI = { 9, "convertMovieToBinaryArray", "/Users/pwjones/Documents/urbanlab/matlab/mousetracking/convertMovieToBinaryArray.m" };
static emlrtRSInfo b_emlrtRSI = { 12, "convertMovieToBinaryArray", "/Users/pwjones/Documents/urbanlab/matlab/mousetracking/convertMovieToBinaryArray.m" };
static emlrtRSInfo c_emlrtRSI = { 13, "convertMovieToBinaryArray", "/Users/pwjones/Documents/urbanlab/matlab/mousetracking/convertMovieToBinaryArray.m" };
static emlrtRSInfo d_emlrtRSI = { 15, "convertMovieToBinaryArray", "/Users/pwjones/Documents/urbanlab/matlab/mousetracking/convertMovieToBinaryArray.m" };
static emlrtRSInfo e_emlrtRSI = { 20, "convertMovieToBinaryArray", "/Users/pwjones/Documents/urbanlab/matlab/mousetracking/convertMovieToBinaryArray.m" };
static emlrtRSInfo f_emlrtRSI = { 23, "convertMovieToBinaryArray", "/Users/pwjones/Documents/urbanlab/matlab/mousetracking/convertMovieToBinaryArray.m" };
static emlrtMCInfo emlrtMCI = { 9, 24, "convertMovieToBinaryArray", "/Users/pwjones/Documents/urbanlab/matlab/mousetracking/convertMovieToBinaryArray.m" };
static emlrtMCInfo b_emlrtMCI = { 15, 16, "convertMovieToBinaryArray", "/Users/pwjones/Documents/urbanlab/matlab/mousetracking/convertMovieToBinaryArray.m" };
static emlrtMCInfo c_emlrtMCI = { 20, 10, "convertMovieToBinaryArray", "/Users/pwjones/Documents/urbanlab/matlab/mousetracking/convertMovieToBinaryArray.m" };
static emlrtMCInfo d_emlrtMCI = { 23, 24, "convertMovieToBinaryArray", "/Users/pwjones/Documents/urbanlab/matlab/mousetracking/convertMovieToBinaryArray.m" };

/* Function Declarations */
static void b_emlrt_marshallIn(const mxArray *u, const emlrtMsgIdentifier *parentId, uint8_T y[307200]);
static void c_emlrt_marshallIn(c_convertMovieToBinaryArrayStac *SD, const mxArray *b_imabsdiff, const char_T *identifier, real_T y[275558400]);
static const mxArray *c_emlrt_marshallOut(const real_T u[275558400]);
static void d_emlrt_marshallIn(c_convertMovieToBinaryArrayStac *SD, const mxArray *u, const emlrtMsgIdentifier *parentId, real_T y[275558400]);
static const mxArray *d_emlrt_marshallOut(const real_T u[307200]);
static void emlrt_marshallIn(const mxArray *b_rgb2gray, const char_T *identifier, uint8_T y[307200]);
static const mxArray *emlrt_marshallOut(const uint8_T u[921600]);
static const mxArray *graythresh(const mxArray *b, emlrtMCInfo *location);
static const mxArray *im2bw(const mxArray *b, const mxArray *c, emlrtMCInfo *location);
static const mxArray *imabsdiff(const mxArray *b, const mxArray *c, emlrtMCInfo *location);
static void m_emlrt_marshallIn(const mxArray *src, const emlrtMsgIdentifier *msgId, uint8_T ret[307200]);
static void n_emlrt_marshallIn(c_convertMovieToBinaryArrayStac *SD, const mxArray *src, const emlrtMsgIdentifier *msgId, real_T ret[275558400]);
static const mxArray *rgb2gray(const mxArray *b, emlrtMCInfo *location);

/* Function Definitions */

static void b_emlrt_marshallIn(const mxArray *u, const emlrtMsgIdentifier *parentId, uint8_T y[307200])
{
    m_emlrt_marshallIn(emlrtAlias(u), parentId, y);
    emlrtDestroyArray(&u);
}

static void c_emlrt_marshallIn(c_convertMovieToBinaryArrayStac *SD, const mxArray *b_imabsdiff, const char_T *identifier, real_T y[275558400])
{
    SD->f1.thisId.fIdentifier = identifier;
    SD->f1.thisId.fParent = NULL;
    d_emlrt_marshallIn(SD, emlrtAlias(b_imabsdiff), &SD->f1.thisId, y);
    emlrtDestroyArray(&b_imabsdiff);
}

static const mxArray *c_emlrt_marshallOut(const real_T u[275558400])
{
    const mxArray *y;
    static const int32_T iv2[3] = { 480, 640, 897 };
    const mxArray *m2;
    real_T (*pData)[];
    int32_T i;
    y = NULL;
    m2 = mxCreateNumericArray(3, (int32_T *)&iv2, mxDOUBLE_CLASS, mxREAL);
    pData = (real_T (*)[])mxGetPr(m2);
    for (i = 0; i < 275558400; i++) {
        (*pData)[i] = u[i];
    }
    emlrtAssign(&y, m2);
    return y;
}

static void d_emlrt_marshallIn(c_convertMovieToBinaryArrayStac *SD, const mxArray *u, const emlrtMsgIdentifier *parentId, real_T y[275558400])
{
    n_emlrt_marshallIn(SD, emlrtAlias(u), parentId, y);
    emlrtDestroyArray(&u);
}

static const mxArray *d_emlrt_marshallOut(const real_T u[307200])
{
    const mxArray *y;
    static const int32_T iv3[2] = { 480, 640 };
    const mxArray *m3;
    real_T (*pData)[];
    int32_T i;
    y = NULL;
    m3 = mxCreateNumericArray(2, (int32_T *)&iv3, mxDOUBLE_CLASS, mxREAL);
    pData = (real_T (*)[])mxGetPr(m3);
    for (i = 0; i < 307200; i++) {
        (*pData)[i] = u[i];
    }
    emlrtAssign(&y, m3);
    return y;
}

static void emlrt_marshallIn(const mxArray *b_rgb2gray, const char_T *identifier, uint8_T y[307200])
{
    emlrtMsgIdentifier thisId;
    thisId.fIdentifier = identifier;
    thisId.fParent = NULL;
    b_emlrt_marshallIn(emlrtAlias(b_rgb2gray), &thisId, y);
    emlrtDestroyArray(&b_rgb2gray);
}

static const mxArray *emlrt_marshallOut(const uint8_T u[921600])
{
    const mxArray *y;
    static const int32_T iv0[3] = { 480, 640, 3 };
    const mxArray *m0;
    uint8_T (*pData)[];
    int32_T i;
    y = NULL;
    m0 = mxCreateNumericArray(3, (int32_T *)&iv0, mxUINT8_CLASS, mxREAL);
    pData = (uint8_T (*)[])mxGetData(m0);
    for (i = 0; i < 921600; i++) {
        (*pData)[i] = u[i];
    }
    emlrtAssign(&y, m0);
    return y;
}

static const mxArray *graythresh(const mxArray *b, emlrtMCInfo *location)
{
    const mxArray *pArray;
    const mxArray *m6;
    pArray = b;
    return emlrtCallMATLAB(1, &m6, 1, &pArray, "graythresh", TRUE, location);
}

static const mxArray *im2bw(const mxArray *b, const mxArray *c, emlrtMCInfo *location)
{
    const mxArray *pArrays[2];
    const mxArray *m7;
    pArrays[0] = b;
    pArrays[1] = c;
    return emlrtCallMATLAB(1, &m7, 2, pArrays, "im2bw", TRUE, location);
}

static const mxArray *imabsdiff(const mxArray *b, const mxArray *c, emlrtMCInfo *location)
{
    const mxArray *pArrays[2];
    const mxArray *m5;
    pArrays[0] = b;
    pArrays[1] = c;
    return emlrtCallMATLAB(1, &m5, 2, pArrays, "imabsdiff", TRUE, location);
}

static void m_emlrt_marshallIn(const mxArray *src, const emlrtMsgIdentifier *msgId, uint8_T ret[307200])
{
    int32_T iv5[2];
    int32_T i3;
    int32_T i4;
    for (i3 = 0; i3 < 2; i3++) {
        iv5[i3] = 480 + 160 * i3;
    }
    emlrtCheckBuiltInCtxR2011b(&emlrtContextGlobal, msgId, src, "uint8", FALSE, 2U, iv5);
    for (i3 = 0; i3 < 640; i3++) {
        for (i4 = 0; i4 < 480; i4++) {
            ret[i4 + 480 * i3] = (*(uint8_T (*)[307200])mxGetData(src))[i4 + 480 * i3];
        }
    }
    emlrtDestroyArray(&src);
}

static void n_emlrt_marshallIn(c_convertMovieToBinaryArrayStac *SD, const mxArray *src, const emlrtMsgIdentifier *msgId, real_T ret[275558400])
{
    int32_T i;
    static const int16_T iv7[3] = { 480, 640, 897 };
    int32_T i5;
    int32_T i6;
    for (i = 0; i < 3; i++) {
        SD->f0.iv6[i] = iv7[i];
    }
    emlrtCheckBuiltInCtxR2011b(&emlrtContextGlobal, msgId, src, "double", FALSE, 3U, SD->f0.iv6);
    for (i = 0; i < 897; i++) {
        for (i5 = 0; i5 < 640; i5++) {
            for (i6 = 0; i6 < 480; i6++) {
                ret[(i6 + 480 * i5) + 307200 * i] = (*(real_T (*)[275558400])mxGetData(src))[(i6 + 480 * i5) + 307200 * i];
            }
        }
    }
    emlrtDestroyArray(&src);
}

static const mxArray *rgb2gray(const mxArray *b, emlrtMCInfo *location)
{
    const mxArray *pArray;
    const mxArray *m4;
    pArray = b;
    return emlrtCallMATLAB(1, &m4, 1, &pArray, "rgb2gray", TRUE, location);
}

void convertMovieToBinaryArray(c_convertMovieToBinaryArrayStac *SD, const b_struct_T mov_struct, uint8_T binArray[275558400])
{
    int32_T ii;
    real_T d0;
    uint8_T u0;
    const mxArray *thresh;
    int16_T b_ii;
    int32_T i0;
    static const mxArray * const threshInitVal = NULL;
    thresh = threshInitVal;
    /*  declare rgb2gray extrinsic, meaning that it stays a matlab function */
    memset((void *)&SD->f2.new_mov[0], 0, 275558400U * sizeof(uint8_T));
    for (ii = 0; ii < 897; ii++) {
        EMLRTPUSHRTSTACK(&emlrtRSI);
        emlrt_marshallIn(rgb2gray(emlrt_marshallOut(mov_struct.frames[ii].cdata), &emlrtMCI), "rgb2gray", *(uint8_T (*)[307200])&SD->f2.new_mov[307200 * ii]);
        EMLRTPOPRTSTACK(&emlrtRSI);
        emlrtBreakCheck();
    }
    /*  make a movie from the average frame to subtract */
    EMLRTPUSHRTSTACK(&b_emlrtRSI);
    mean(SD->f2.new_mov, SD->f2.dv0);
    b_round(SD->f2.dv0);
    for (ii = 0; ii < 307200; ii++) {
        d0 = SD->f2.dv0[ii];
        d0 = d0 < 0.0 ? muDoubleScalarCeil(d0 - 0.5) : muDoubleScalarFloor(d0 + 0.5);
        if (d0 < 256.0) {
            if (d0 >= 0.0) {
                u0 = (uint8_T)d0;
            } else {
                u0 = 0;
            }
        } else if (d0 >= 256.0) {
            u0 = MAX_uint8_T;
        } else {
            u0 = 0;
        }
        SD->f2.avg_frame[ii] = u0;
    }
    EMLRTPOPRTSTACK(&b_emlrtRSI);
    EMLRTPUSHRTSTACK(&c_emlrtRSI);
    repmat(SD->f2.avg_frame, SD->f2.avg_mov);
    EMLRTPOPRTSTACK(&c_emlrtRSI);
    /*  need to declare the array so that the MATLAB lets me index it later */
    EMLRTPUSHRTSTACK(&d_emlrtRSI);
    c_emlrt_marshallIn(SD, imabsdiff(b_emlrt_marshallOut(SD->f2.new_mov), b_emlrt_marshallOut(SD->f2.avg_mov), &b_emlrtMCI), "imabsdiff", SD->f2.mean_sub_mov);
    EMLRTPOPRTSTACK(&d_emlrtRSI);
    /* this should give a nice moving blob. */
    memset((void *)&binArray[0], 0, 275558400U * sizeof(uint8_T));
    EMLRTPUSHRTSTACK(&e_emlrtRSI);
    emlrtAssign(&thresh, graythresh(c_emlrt_marshallOut(SD->f2.mean_sub_mov), &c_emlrtMCI));
    EMLRTPOPRTSTACK(&e_emlrtRSI);
    for (b_ii = 0; b_ii < 897; b_ii++) {
        EMLRTPUSHRTSTACK(&f_emlrtRSI);
        emlrt_marshallIn(im2bw(d_emlrt_marshallOut(*(real_T (*)[307200])&SD->f2.mean_sub_mov[307200 * b_ii]), emlrtAlias(thresh), &d_emlrtMCI), "im2bw", *(uint8_T (*)[307200])&binArray[307200 * b_ii]);
        EMLRTPOPRTSTACK(&f_emlrtRSI);
        emlrtBreakCheck();
    }
    for (ii = 0; ii < 275558400; ii++) {
        i0 = (int32_T)((uint32_T)binArray[ii] * 255U);
        if ((uint32_T)i0 > 255U) {
            i0 = 255;
        }
        binArray[ii] = (uint8_T)i0;
    }
    emlrtDestroyArray(&thresh);
}
/* End of code generation (convertMovieToBinaryArray.c) */
