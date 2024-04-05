#include <stdio.h>

#define FRAME_SIZE 8
#define FILT_LEN 3

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

int main()
{
    float x[32] = { 1, 0, 0, 0,0, 1, 0, 0, 0,0, 1, 0,0,0,0,0,0,0,1,1,1,1,1,1,0,0,0,0,0,0,0,0 }, * y_ptr, * x_ptr;
    float h[FILT_LEN] = { 1, 2, 3 };
    float y[40] = { 0 }, temp, out[40] = { 0 };
    int i = 0, j,len;
    len = (FRAME_SIZE + FILT_LEN - 1);
    for (i = 0; i < 4; i++)
    {
        y_ptr = &y[i * len];
        x_ptr = &x[i * FRAME_SIZE];
        fir_filter(x_ptr, h, y_ptr, FRAME_SIZE, FILT_LEN);
    }

    for(i = 0; i < FRAME_SIZE; i++)
    {
        out[i] = y[i];
    }

    y_ptr = &y[0];
    for (j = 1; j < 4; j++)
    {
        for (i = 0; i < (FILT_LEN-1); i++)
        {
            out[j* FRAME_SIZE +i] = y_ptr[len*(j-1) + FRAME_SIZE + i] + y_ptr[len * (j - 1) + len + i];
        }
        for (i = 0; i < (FRAME_SIZE - FILT_LEN + 1); i++)
        {
            out[j * FRAME_SIZE + (FILT_LEN - 1 + i)] = y_ptr[(j)*len + (FILT_LEN - 1 + i)];
        }
    }
    for (j = 0; j < 26; j++)
        printf("%f,%f\n", out[j],y[j]);
    return 0;
}