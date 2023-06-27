#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define USE_FIXED_PT_CODE

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
    return (((float)x)/(pow(2,qfactor)));
}

void fir_filter(float* in, float* coeffs, float* out, float *zi,int num_of_filt_coeffs, int frame_size)
{
	float delay_line[511];
	int i,j,count_r,count_f = 0,index;
	float sum;

	for (i = 0; i < (num_of_filt_coeffs-1); i++)
	{
		delay_line[i+1] = zi[i];
	}

	count_f = 0;
	for (i = 0; i < frame_size; i++)
	{
		if (count_f >= num_of_filt_coeffs)
		{
			index = 0;
			count_f = 0;
		}
		else
		{
			index = count_f;
		}
		delay_line[index] = in[i];
		count_r = index;
		sum = 0;
		for (j = 0; j < num_of_filt_coeffs; j++)
		{
			if (count_r < 0)
			{
				index = count_r + num_of_filt_coeffs;
			}
			else
			{
				index = count_r;
			}
			sum += coeffs[j] * delay_line[index];
			count_r--;
		}
		out[i] = sum;
		count_f++;
	}
	index = frame_size - num_of_filt_coeffs + 1;
	for (i = 0; i < (num_of_filt_coeffs-1); i++)
	{
		zi[i] = in[index];
		index++;
	}
}

void fir_filter_fxd_pt(Word32* in, Word32* coeffs, Word32* out, Word32 *zi,Word32 num_of_filt_coeffs, Word32 frame_size)
{
	Word32 delay_line[511];
	Word32 i,j,count_r,count_f = 0,index;
	Word64 sum;

	for (i = 0; i < (num_of_filt_coeffs-1); i++)
	{
		delay_line[i+1] = zi[i];
	}

	count_f = 0;
	for (i = 0; i < frame_size; i++)
	{
		if (count_f >= num_of_filt_coeffs)
		{
			index = 0;
			count_f = 0;
		}
		else
		{
			index = count_f;
		}
		delay_line[index] = in[i];
		count_r = index;
		sum = 0;
		for (j = 0; j < num_of_filt_coeffs; j++)
		{
			if (count_r < 0)
			{
				index = count_r + num_of_filt_coeffs;
			}
			else
			{
				index = count_r;
			}
			//sum += coeffs[j] * delay_line[index];
			sum = s64_mla_s32_s32(sum, coeffs[j], delay_line[index]); //Q2.29*Q2.29 = Q4.58
			count_r--;
		}
		out[i] = (Word32)(sum >> 31); //Q4.27
		count_f++;
	}
	index = frame_size - num_of_filt_coeffs + 1;
	for (i = 0; i < (num_of_filt_coeffs-1); i++)
	{
		zi[i] = in[index];
		index++;
	}
}

int main(void)
{
	FILE *fcoeffs, *finput, *fout;
	float in[4000], coeffs[511], zi[510], out[4000];
	Word32 coeffs_fxd_pt[511],out_fxd_pt[4000];
	Word32 in_fxd_pt[4000],zi_fxd_pt[510];
	int i,j;

	fcoeffs = fopen("fir_ceoffs_pygen.bin","rb");
	finput = fopen("input_pygen.bin", "rb");
	fout = fopen("out_msvc.bin","wb");

	fread(coeffs,511,sizeof(float),fcoeffs);
	for (i = 0; i < 511; i++)
    {
        coeffs_fxd_pt[i] = float_to_fixed_conv(coeffs[i],29);
    }

	for (i = 0; i < 510; i++)
	{
		zi[i] = 0.0f;
		zi_fxd_pt[i] = 0;
	}

	for (j = 0; j < 4; j++)
	{
		fread(in, 4000, sizeof(float), finput);
		for (i = 0; i < 4000; i++)
        {
           in_fxd_pt[i] = float_to_fixed_conv(in[i],29);
        }
#ifdef USE_FIXED_PT_CODE
		fir_filter_fxd_pt(in_fxd_pt, coeffs_fxd_pt, out_fxd_pt, zi_fxd_pt, 511, 4000);
		for (i = 0; i < 4000; i++)
        {
           out[i] = fixed_to_float_conv(out_fxd_pt[i],27);
        }
#else
        fir_filter(in, coeffs, out, zi, 511, 4000);
#endif
		fwrite(out,4000,sizeof(float),fout);
	}

	fclose(fcoeffs);
	fclose(finput);
	fclose(fout);

	return 0;
}
