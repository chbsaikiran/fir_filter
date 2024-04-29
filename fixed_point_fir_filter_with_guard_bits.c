#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define ABS_FLOAT(x) ((x) > 0 ? (x):(0-(x)))

short float_to_fixed_conv_16bit(float x, int qfactor)
{
    return ((short)(x*(pow(2,qfactor))));
}

float fixed_to_float_conv_16bit(short x, int qfactor)
{
    return (((float)x)/(pow(2,qfactor)));
}

void fir_filter_fixed(short *x,short *h,short *y,int N,int M)
{
    int n,k;
    long long sum;
    for(n = 0; n < (N+M-1); n++)
    {
        sum = 0;
        for(k = 0; k < M; k++)
        {
            if((n-k) >=0 && (n-k) < N)
            {
                sum += ((int)x[n-k])*((int)h[k]);
            }
        }
        y[n] = (short)(sum >> 17); //Q4.11 //17 is got from (15 + 2), 2 is got from log2(0.5*8)
    }
}

void fir_filter(float* x, float* h, float* y, int N, int M)
{
    int n, k;
    float sum;
    for (n = 0; n < (N + M - 1); n++)
    {
        sum = 0.0f;
        for (k = 0; k < M; k++)
        {
            if ((n - k) >= 0 && (n - k) < N)
            {
                sum += x[n - k] * h[k];
            }
        }
        y[n] = sum;
    }
}

/* 
   1. The output Q format will be Q4.11
   2. If Guard bit were not there then the Q format of h would be Q2.13
   3. Because of Guard bits we were able to use a Q format of Q0.15 for h 
   4. We have used more precision for h due to guard bits, that is what this code 
      is trying to establish
 */

int main()
{
    float h_float[8] = { 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f };
    float in_float[16], y_float[23],y_fix_to_float[23], max_error = -1.0f;
    int i, qf_h = 15, qf_in = 13,qf_out = (15 - (0+(15-qf_in)+2)); //2 is got from log2(0.5*8)
    short in[16],h[8],y[23];

    for (i = 0; i < 16; i += 2)
    {
        in_float[i] = 1.75f;
        in_float[i + 1] = 1.5f;
    }

    for (i = 0; i < 16; i++)
    {
        in[i] = float_to_fixed_conv_16bit(in_float[i], qf_in);
    }

    for (i = 0; i < 8; i++)
    {
        h[i] = float_to_fixed_conv_16bit(h_float[i], qf_h);
    }
    fir_filter_fixed(in,h,y,16,8);
    fir_filter(in_float, h_float, y_float, 16, 8);

    for (i = 0; i < 23; i++)
    {
        y_fix_to_float[i] = fixed_to_float_conv_16bit(y[i], qf_out);

        printf("y_fix_to_float[%d] = %f,y_float[%d] = %f\n", i,y_fix_to_float[i],i, y_float[i]);
        if (ABS_FLOAT((y_float[i]-y_fix_to_float[i])) > max_error)
        {
            max_error = ABS_FLOAT((y_float[i] - y_fix_to_float[i]));
        }
    }

    printf("max_error = %f\n", max_error);

    return 0;
}