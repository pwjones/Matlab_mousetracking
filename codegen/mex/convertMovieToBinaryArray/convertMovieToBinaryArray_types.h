/*
 * convertMovieToBinaryArray_types.h
 *
 * Code generation for function 'convertMovieToBinaryArray'
 *
 * C source code generated on: Tue Apr 24 12:18:22 2012
 *
 */

#ifndef __CONVERTMOVIETOBINARYARRAY_TYPES_H__
#define __CONVERTMOVIETOBINARYARRAY_TYPES_H__

/* Type Definitions */
typedef struct
{
    uint8_T cdata[921600];
} struct_T;
typedef struct
{
    real_T width;
    real_T height;
    real_T nrFramesTotal;
    struct_T frames[897];
    real_T rate;
    real_T totalDuration;
    real_T times[897];
    boolean_T skippedFrames;
} b_struct_T;
typedef struct
{
    struct
    {
        int32_T iv6[3];
    } f0;
    struct
    {
        emlrtMsgIdentifier thisId;
    } f1;
    struct
    {
        uint8_T new_mov[275558400];
        uint8_T avg_mov[275558400];
        real_T dv0[307200];
        uint8_T avg_frame[307200];
        real_T mean_sub_mov[275558400];
    } f2;
    struct
    {
        b_struct_T expl_temp;
        struct_T mov_struct_frames[897];
        uint8_T binArray[275558400];
        real_T mov_struct_times[897];
    } f3;
} c_convertMovieToBinaryArrayStac;

#endif
/* End of code generation (convertMovieToBinaryArray_types.h) */
