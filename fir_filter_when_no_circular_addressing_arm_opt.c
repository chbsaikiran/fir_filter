#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "arm_neon.h"

#define USE_FIXED_PT_CODE
#define PROFILE_CODE

#ifdef PROFILE_CODE
#include <sys/time.h>
#endif

#define ABS_FLOAT(x) ((x) > 0 ? (x):(0-(x)))

#define Word64 long long
#define Word32 int
#define Word16 short

Word64 s64_mul_s32_s32(Word32 x, Word32 y)
{
    Word64 prod;

    prod = ((Word64)x)*((Word64)y);

    return prod;
}

Word64 s64_mla_s32_s32(Word64 sum,Word32 x, Word32 y)
{
    Word64 prod;

    prod = ((Word64)x)*((Word64)y);

    sum = sum + prod;

    return sum;
}

Word32 float_to_fixed_conv(float x, Word32 qfactor)
{
    return ((Word32)(x*(pow(2,qfactor))));
}

Word16 float_to_fixed_conv_16bit(float x, Word32 qfactor)
{
    return ((Word16)(x*(pow(2,qfactor))));
}

float fixed_to_float_conv(Word32 x, Word32 qfactor)
{
    return (((float)x)/((float)(pow(2,qfactor))));
}

Word32 delay_line[510+4000];

void fir_filter_fxd_pt(Word32* in, Word32* coeffs, Word32* out,Word32 num_of_filt_coeffs, Word32 frame_size)
{
    Word64 acc1,acc2;     // accumulator for MACs
    Word32 *coeffp; // pointer to coefficients
    Word32 *inputp,*inputp1; // pointer to input samples
    Word32 n,len,rem;
    Word32 k;
    int32x4_t q0,q1,q2,q3;
    int64x2_t q5,q4,q6,q7;
    int64x1_t d12,d13,d14,d15;
    int32x2_t d0;
 
    // put the new samples at the high end of the buffer
    memcpy( &delay_line[num_of_filt_coeffs - 1], in,
            frame_size * sizeof(int) );
 
    // apply the filter to each input sample
    for ( n = 0; n < frame_size; n+=2 ) 
    {
        // calculate out n
        coeffp = &coeffs[num_of_filt_coeffs-1];
        inputp = &delay_line[num_of_filt_coeffs + n];
        inputp1 = inputp - 1;
        acc1 = 0;
        acc2 = 0;
        len = num_of_filt_coeffs >> 2;
        rem = num_of_filt_coeffs & 3;
        while(rem--)
        {
            acc1 = s64_mla_s32_s32(acc1, (*coeffp), (*inputp--));
            acc2 = s64_mla_s32_s32(acc2, (*coeffp), (*inputp1--));
            coeffp--;
        }
        d14 = vdup_n_s64(acc2);
        d15 = vdup_n_s64(acc1);
        q7 = vcombine_s64(d14,d15);
        coeffp = coeffp - 3;
        inputp = inputp - 3;
        inputp1 = inputp1 - 3;
        q4 = vdupq_n_s64(0);
        q5 = vdupq_n_s64(0);
        while(len--) 
        {
            //acc = s64_mla_s32_s32(acc, (*coeffp--), (*inputp--));
            q0 = vld1q_s32(((int32_t*)coeffp)); //3 2 1 0
            q1 = vld1q_s32(((int32_t*)inputp)); //7 6 5 4
            inputp -= 4;
            coeffp -= 4;
            q2 = vld1q_s32(((int32_t*)inputp1)); //6 5 4 3
            inputp1 -= 4;
            q5 = vmlal_s32(q5,vget_low_s32(q0),vget_low_s32(q1));
            q4 = vmlal_s32(q4,vget_low_s32(q0),vget_low_s32(q2));
            q5 = vmlal_s32(q5,vget_high_s32(q0),vget_high_s32(q1));
            q4 = vmlal_s32(q4,vget_high_s32(q0),vget_high_s32(q2));
        }
        d12 = vadd_s64(vget_low_s64(q4),vget_high_s64(q4));
        d13 = vadd_s64(vget_low_s64(q5),vget_high_s64(q5));
        q6 = vcombine_s64(d12,d13);
        q6 = vaddq_s64(q7,q6);
        d0 = vshrn_n_s64(q6,31);
        //out[n] = (Word32)(acc >> 31);
        vst1_s32(((int32_t*)&out[n]), d0);
    }
    // shift input samples back in time for next time
    memcpy( &delay_line[0], &delay_line[frame_size],
            (num_of_filt_coeffs - 1) * sizeof(int));
 
}

int main(void)
{
    FILE *fcoeffs, *finput, *fout;
    float in[4000], coeffs[511],out[4000];
    Word32 coeffs_fxd_pt[511],out_fxd_pt[4000];
    Word32 in_fxd_pt[4000],temp;
    int i,j;
#ifdef PROFILE_CODE
    long seconds;
    long microseconds;
    double elapsed = 0;
#endif

    fcoeffs = fopen("fir_ceoffs_pygen.bin","rb");
    finput = fopen("input_pygen.bin", "rb");
    fout = fopen("out_arm_wo_circ_buffer.bin","wb");

    fread(coeffs,511,sizeof(float),fcoeffs);
    for (i = 0; i < 511; i++)
    {
        coeffs_fxd_pt[i] = float_to_fixed_conv(coeffs[i],29);
    }
    for (i = 0; i < 511; i++)
    {
        temp = coeffs_fxd_pt[i];
        coeffs_fxd_pt[i] = coeffs_fxd_pt[511 - i];
        coeffs_fxd_pt[511 - i] = temp;
    }

    for (i = 0; i < (510+4000); i++)
    {
        delay_line[i] = 0;
    }
    
    for (j = 0; j < 4; j++)
    {
        fread(in, 4000, sizeof(float), finput);
        for (i = 0; i < 4000; i++)
        {
           in_fxd_pt[i] = float_to_fixed_conv(in[i],29);
        }
#ifdef PROFILE_CODE
        struct timeval start, end;
        gettimeofday(&start, NULL);
#endif
        fir_filter_fxd_pt(in_fxd_pt, coeffs_fxd_pt, out_fxd_pt, 511, 4000);
#ifdef PROFILE_CODE
        gettimeofday(&end, NULL);
        seconds = (end.tv_sec - start.tv_sec);
        microseconds = ((seconds * 1000000) + end.tv_usec) - (start.tv_usec);
        elapsed += microseconds*1e-6;
#endif
        for (i = 0; i < 4000; i++)
        {
           out[i] = fixed_to_float_conv(out_fxd_pt[i],27);
        }
        fwrite(out,4000,sizeof(float),fout);
    }
#ifdef PROFILE_CODE
    printf("elapsed_time = %lf\n",elapsed);
#endif

    fclose(fcoeffs);
    fclose(finput);
    fclose(fout);

    return 0;
}
