/*
 * convertMovieToBinaryArray_api.c
 *
 * Code generation for function 'convertMovieToBinaryArray_api'
 *
 * C source code generated on: Tue Apr 24 12:18:22 2012
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "convertMovieToBinaryArray.h"
#include "convertMovieToBinaryArray_api.h"
#include "convertMovieToBinaryArray_mexutil.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */
static void e_emlrt_marshallIn(const mxArray *mov_struct, const char_T *identifier, real_T *y_width, real_T *y_height, real_T *y_nrFramesTotal, struct_T y_frames[897], real_T *y_rate, real_T *y_totalDuration, real_T y_times[897], boolean_T *y_skippedFrames);
static void f_emlrt_marshallIn(const mxArray *u, const emlrtMsgIdentifier *parentId, real_T *y_width, real_T *y_height, real_T *y_nrFramesTotal, struct_T y_frames[897], real_T *y_rate, real_T *y_totalDuration, real_T y_times[897], boolean_T *y_skippedFrames);
static real_T g_emlrt_marshallIn(const mxArray *u, const emlrtMsgIdentifier *parentId);
static void h_emlrt_marshallIn(const mxArray *u, const emlrtMsgIdentifier *parentId, struct_T y[897]);
static void i_emlrt_marshallIn(const mxArray *u, const emlrtMsgIdentifier *parentId, uint8_T y[921600]);
static void j_emlrt_marshallIn(const mxArray *u, const emlrtMsgIdentifier *parentId);
static void k_emlrt_marshallIn(const mxArray *u, const emlrtMsgIdentifier *parentId, real_T y[897]);
static boolean_T l_emlrt_marshallIn(const mxArray *u, const emlrtMsgIdentifier *parentId);
static real_T o_emlrt_marshallIn(const mxArray *src, const emlrtMsgIdentifier *msgId);
static void p_emlrt_marshallIn(const mxArray *src, const emlrtMsgIdentifier *msgId, uint8_T ret[921600]);
static void q_emlrt_marshallIn(const mxArray *src, const emlrtMsgIdentifier *msgId);
static void r_emlrt_marshallIn(const mxArray *src, const emlrtMsgIdentifier *msgId, real_T ret[897]);
static boolean_T s_emlrt_marshallIn(const mxArray *src, const emlrtMsgIdentifier *msgId);

/* Function Definitions */

static void e_emlrt_marshallIn(const mxArray *mov_struct, const char_T *identifier, real_T *y_width, real_T *y_height, real_T *y_nrFramesTotal, struct_T y_frames[897], real_T *y_rate, real_T *y_totalDuration, real_T y_times[897], boolean_T *y_skippedFrames)
{
    emlrtMsgIdentifier thisId;
    thisId.fIdentifier = identifier;
    thisId.fParent = NULL;
    f_emlrt_marshallIn(emlrtAlias(mov_struct), &thisId, y_width, y_height, y_nrFramesTotal, y_frames, y_rate, y_totalDuration, y_times, y_skippedFrames);
    emlrtDestroyArray(&mov_struct);
}

static void f_emlrt_marshallIn(const mxArray *u, const emlrtMsgIdentifier *parentId, real_T *y_width, real_T *y_height, real_T *y_nrFramesTotal, struct_T y_frames[897], real_T *y_rate, real_T *y_totalDuration, real_T y_times[897], boolean_T *y_skippedFrames)
{
    emlrtMsgIdentifier thisId;
    static const char * fieldNames[8] = { "width", "height", "nrFramesTotal", "frames", "rate", "totalDuration", "times", "skippedFrames" };
    thisId.fParent = parentId;
    emlrtCheckStructR2011a(parentId, u, 8, fieldNames, 0U, 0);
    thisId.fIdentifier = "width";
    *y_width = g_emlrt_marshallIn(emlrtAlias(emlrtGetField(u, 0, "width")), &thisId);
    thisId.fIdentifier = "height";
    *y_height = g_emlrt_marshallIn(emlrtAlias(emlrtGetField(u, 0, "height")), &thisId);
    thisId.fIdentifier = "nrFramesTotal";
    *y_nrFramesTotal = g_emlrt_marshallIn(emlrtAlias(emlrtGetField(u, 0, "nrFramesTotal")), &thisId);
    thisId.fIdentifier = "frames";
    h_emlrt_marshallIn(emlrtAlias(emlrtGetField(u, 0, "frames")), &thisId, y_frames);
    thisId.fIdentifier = "rate";
    *y_rate = g_emlrt_marshallIn(emlrtAlias(emlrtGetField(u, 0, "rate")), &thisId);
    thisId.fIdentifier = "totalDuration";
    *y_totalDuration = g_emlrt_marshallIn(emlrtAlias(emlrtGetField(u, 0, "totalDuration")), &thisId);
    thisId.fIdentifier = "times";
    k_emlrt_marshallIn(emlrtAlias(emlrtGetField(u, 0, "times")), &thisId, y_times);
    thisId.fIdentifier = "skippedFrames";
    *y_skippedFrames = l_emlrt_marshallIn(emlrtAlias(emlrtGetField(u, 0, "skippedFrames")), &thisId);
    emlrtDestroyArray(&u);
}

static real_T g_emlrt_marshallIn(const mxArray *u, const emlrtMsgIdentifier *parentId)
{
    real_T y;
    y = o_emlrt_marshallIn(emlrtAlias(u), parentId);
    emlrtDestroyArray(&u);
    return y;
}

static void h_emlrt_marshallIn(const mxArray *u, const emlrtMsgIdentifier *parentId, struct_T y[897])
{
    emlrtMsgIdentifier thisId;
    int32_T iv4[2];
    int32_T i2;
    static const char * fieldNames[2] = { "cdata", "colormap" };
    struct_T (*r0)[897];
    thisId.fParent = parentId;
    for (i2 = 0; i2 < 2; i2++) {
        iv4[i2] = 1 + 896 * i2;
    }
    emlrtCheckStructR2011a(parentId, u, 2, fieldNames, 2U, iv4);
    r0 = (struct_T (*)[897])y;
    for (i2 = 0; i2 < 897; i2++) {
        thisId.fIdentifier = "cdata";
        i_emlrt_marshallIn(emlrtAlias(emlrtGetField(u, i2, "cdata")), &thisId, (*r0)[i2].cdata);
        thisId.fIdentifier = "colormap";
        j_emlrt_marshallIn(emlrtAlias(emlrtGetField(u, i2, "colormap")), &thisId);
    }
    emlrtDestroyArray(&u);
}

static void i_emlrt_marshallIn(const mxArray *u, const emlrtMsgIdentifier *parentId, uint8_T y[921600])
{
    p_emlrt_marshallIn(emlrtAlias(u), parentId, y);
    emlrtDestroyArray(&u);
}

static void j_emlrt_marshallIn(const mxArray *u, const emlrtMsgIdentifier *parentId)
{
    q_emlrt_marshallIn(emlrtAlias(u), parentId);
    emlrtDestroyArray(&u);
}

static void k_emlrt_marshallIn(const mxArray *u, const emlrtMsgIdentifier *parentId, real_T y[897])
{
    r_emlrt_marshallIn(emlrtAlias(u), parentId, y);
    emlrtDestroyArray(&u);
}

static boolean_T l_emlrt_marshallIn(const mxArray *u, const emlrtMsgIdentifier *parentId)
{
    boolean_T y;
    y = s_emlrt_marshallIn(emlrtAlias(u), parentId);
    emlrtDestroyArray(&u);
    return y;
}

static real_T o_emlrt_marshallIn(const mxArray *src, const emlrtMsgIdentifier *msgId)
{
    real_T ret;
    emlrtCheckBuiltInCtxR2011b(&emlrtContextGlobal, msgId, src, "double", FALSE, 0U, 0);
    ret = *(real_T *)mxGetData(src);
    emlrtDestroyArray(&src);
    return ret;
}

static void p_emlrt_marshallIn(const mxArray *src, const emlrtMsgIdentifier *msgId, uint8_T ret[921600])
{
    int32_T iv8[3];
    int32_T i;
    static const int16_T iv9[3] = { 480, 640, 3 };
    int32_T i7;
    int32_T i8;
    for (i = 0; i < 3; i++) {
        iv8[i] = iv9[i];
    }
    emlrtCheckBuiltInCtxR2011b(&emlrtContextGlobal, msgId, src, "uint8", FALSE, 3U, iv8);
    for (i = 0; i < 3; i++) {
        for (i7 = 0; i7 < 640; i7++) {
            for (i8 = 0; i8 < 480; i8++) {
                ret[(i8 + 480 * i7) + 307200 * i] = (*(uint8_T (*)[921600])mxGetData(src))[(i8 + 480 * i7) + 307200 * i];
            }
        }
    }
    emlrtDestroyArray(&src);
}

static void q_emlrt_marshallIn(const mxArray *src, const emlrtMsgIdentifier *msgId)
{
    int32_T iv10[2];
    int32_T i;
    for (i = 0; i < 2; i++) {
        iv10[i] = 0;
    }
    emlrtCheckBuiltInCtxR2011b(&emlrtContextGlobal, msgId, src, "double", FALSE, 2U, iv10);
    emlrtDestroyArray(&src);
}

static void r_emlrt_marshallIn(const mxArray *src, const emlrtMsgIdentifier *msgId, real_T ret[897])
{
    int32_T iv11[2];
    int32_T i9;
    for (i9 = 0; i9 < 2; i9++) {
        iv11[i9] = 1 + 896 * i9;
    }
    emlrtCheckBuiltInCtxR2011b(&emlrtContextGlobal, msgId, src, "double", FALSE, 2U, iv11);
    for (i9 = 0; i9 < 897; i9++) {
        ret[i9] = (*(real_T (*)[897])mxGetData(src))[i9];
    }
    emlrtDestroyArray(&src);
}

static boolean_T s_emlrt_marshallIn(const mxArray *src, const emlrtMsgIdentifier *msgId)
{
    boolean_T ret;
    emlrtCheckBuiltInCtxR2011b(&emlrtContextGlobal, msgId, src, "logical", FALSE, 0U, 0);
    ret = *mxGetLogicals(src);
    emlrtDestroyArray(&src);
    return ret;
}

void convertMovieToBinaryArray_api(c_convertMovieToBinaryArrayStac *SD, const mxArray * const prhs[1], const mxArray *plhs[1])
{
    boolean_T mov_struct_skippedFrames;
    real_T mov_struct_totalDuration;
    real_T mov_struct_rate;
    real_T mov_struct_nrFramesTotal;
    real_T mov_struct_height;
    real_T mov_struct_width;
    int32_T i1;
    /* Marshall function inputs */
    e_emlrt_marshallIn(emlrtAliasP(prhs[0]), "mov_struct", &mov_struct_width, &mov_struct_height, &mov_struct_nrFramesTotal, SD->f3.mov_struct_frames, &mov_struct_rate, &mov_struct_totalDuration, SD->f3.mov_struct_times, &mov_struct_skippedFrames);
    /* Invoke the target function */
    SD->f3.expl_temp.skippedFrames = mov_struct_skippedFrames;
    for (i1 = 0; i1 < 897; i1++) {
        SD->f3.expl_temp.times[i1] = SD->f3.mov_struct_times[i1];
    }
    SD->f3.expl_temp.totalDuration = mov_struct_totalDuration;
    SD->f3.expl_temp.rate = mov_struct_rate;
    for (i1 = 0; i1 < 897; i1++) {
        SD->f3.expl_temp.frames[i1] = SD->f3.mov_struct_frames[i1];
    }
    SD->f3.expl_temp.nrFramesTotal = mov_struct_nrFramesTotal;
    SD->f3.expl_temp.height = mov_struct_height;
    SD->f3.expl_temp.width = mov_struct_width;
    convertMovieToBinaryArray(SD, SD->f3.expl_temp, SD->f3.binArray);
    /* Marshall function outputs */
    plhs[0] = b_emlrt_marshallOut(SD->f3.binArray);
}
/* End of code generation (convertMovieToBinaryArray_api.c) */
