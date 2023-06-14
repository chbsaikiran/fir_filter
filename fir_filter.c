#include <stdio.h>
#include <stdlib.h>

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

int main(void)
{
	FILE *fcoeffs, *finput, *fout;
	float in[4000], coeffs[511], zi[510], out[4000];
	int i;

	fcoeffs = fopen("fir_ceoffs_pygen.bin","rb");
	finput = fopen("input_pygen.bin", "rb");
	fout = fopen("out_msvc.bin","wb");

	fread(coeffs,511,sizeof(float),fcoeffs);
	for (i = 0; i < 510; i++)
	{
		zi[i] = 0.0f;
	}

	for (i = 0; i < 4; i++)
	{
		fread(in, 4000, sizeof(float), finput);
		fir_filter(in, coeffs, out, zi, 511, 4000);
		fwrite(out,4000,sizeof(float),fout);
	}

	fclose(fcoeffs);
	fclose(finput);
	fclose(fout);

	return 0;
}