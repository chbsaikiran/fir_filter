#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define USE_FIXED_PT_CODE
//#define PROFILE_CODE

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
    Word64 acc;     // accumulator for MACs
    Word32 *coeffp; // pointer to coefficients
    Word32 *inputp; // pointer to input samples
    Word32 n;
    Word32 k;
 
    // put the new samples at the high end of the buffer
    memcpy( &delay_line[num_of_filt_coeffs - 1], in,
            frame_size * sizeof(int) );
 
    // apply the filter to each input sample
    for ( n = 0; n < frame_size; n++ ) 
    {
        // calculate out n
        coeffp = coeffs;
        inputp = &delay_line[num_of_filt_coeffs - 1 + n];
        acc = 0;
        for ( k = 0; k < num_of_filt_coeffs; k++ ) 
        {
            acc = s64_mla_s32_s32(acc, (*coeffp++), (*inputp--));
        }
        out[n] = (Word32)(acc >> 31);
    }
    // shift input samples back in time for next time
    memmove( &delay_line[0], &delay_line[frame_size],
            (num_of_filt_coeffs - 1) * sizeof(int));
 
}

int main(void)
{
    FILE *fcoeffs, *finput, *fout;
    float in[4000], coeffs[511],out[4000];
    Word32 coeffs_fxd_pt[511],out_fxd_pt[4000];
    Word32 in_fxd_pt[4000];
    int i,j;
#ifdef PROFILE_CODE
    long seconds;
    long microseconds;
    double elapsed = 0;
#endif

    fcoeffs = fopen("..\\fir_ceoffs_pygen.bin","rb");
    finput = fopen("..\\input_pygen.bin", "rb");
    fout = fopen("..\\out_msvc_wo_circ_buffer.bin","wb");

    fread(coeffs,511,sizeof(float),fcoeffs);
    for (i = 0; i < 511; i++)
    {
        coeffs_fxd_pt[i] = float_to_fixed_conv(coeffs[i],29);
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
