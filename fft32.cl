/* FFT
 * command: python -m fftsynth.generator --radix 2 --depth 5 --fpga --fma --ctype float8
 */
/* ~\~ language=OpenCL filename=fftsynth/templates/preprocessor.cl */
/* ~\~ begin <<lit/code-generator.md|fftsynth/templates/preprocessor.cl>>[0] */
#pragma OPENCL EXTENSION cl_intel_channels : enable
#include <ihc_apint.h>
#ifdef TESTING
channel float8 in_channel, out_channel;
#endif // TESTING
#define SWAP(type, x, y) do { type temp = x; x = y, y = temp; } while ( false );

#define DIVR(x) ((x) >> 1)
#define MODR(x) ((x) & 1)
#define MULR(x) ((x) << 1)

/* ~\~ end */
/* ~\~ language=OpenCL filename=fftsynth/templates/fma-twiddles.cl */
/* ~\~ begin <<lit/fma-codelets.md|fftsynth/templates/fma-twiddles.cl>>[0] */
__constant float2 W[64][1] = {

{ (float2)(1.000000f, 0.000000f) },
{ (float2)(0.000000f, -1.000000f) },
{ (float2)(0.000000f, -1.000000f) },
{ (float2)(1.000000f, 0.000000f) },
{ (float2)(0.000000f, -1.000000f) },
{ (float2)(1.000000f, 0.000000f) },
{ (float2)(1.000000f, 0.000000f) },
{ (float2)(0.000000f, -1.000000f) },
{ (float2)(0.000000f, -1.000000f) },
{ (float2)(1.000000f, 0.000000f) },
{ (float2)(1.000000f, 0.000000f) },
{ (float2)(0.000000f, -1.000000f) },
{ (float2)(1.000000f, 0.000000f) },
{ (float2)(0.000000f, -1.000000f) },
{ (float2)(0.000000f, -1.000000f) },
{ (float2)(1.000000f, 0.000000f) },
{ (float2)(1.000000f, 0.000000f) },
{ (float2)(0.707107f, -0.707107f) },
{ (float2)(-0.707107f, -0.707107f) },
{ (float2)(0.000000f, -1.000000f) },
{ (float2)(0.707107f, -0.707107f) },
{ (float2)(1.000000f, 0.000000f) },
{ (float2)(0.000000f, -1.000000f) },
{ (float2)(-0.707107f, -0.707107f) },
{ (float2)(0.707107f, -0.707107f) },
{ (float2)(1.000000f, 0.000000f) },
{ (float2)(0.000000f, -1.000000f) },
{ (float2)(-0.707107f, -0.707107f) },
{ (float2)(1.000000f, 0.000000f) },
{ (float2)(0.707107f, -0.707107f) },
{ (float2)(-0.707107f, -0.707107f) },
{ (float2)(0.000000f, -1.000000f) },
{ (float2)(1.000000f, 0.000000f) },
{ (float2)(0.923880f, -0.382683f) },
{ (float2)(0.382683f, -0.923880f) },
{ (float2)(0.707107f, -0.707107f) },
{ (float2)(-0.382683f, -0.923880f) },
{ (float2)(0.000000f, -1.000000f) },
{ (float2)(-0.707107f, -0.707107f) },
{ (float2)(-0.923880f, -0.382683f) },
{ (float2)(0.923880f, -0.382683f) },
{ (float2)(1.000000f, 0.000000f) },
{ (float2)(0.707107f, -0.707107f) },
{ (float2)(0.382683f, -0.923880f) },
{ (float2)(0.000000f, -1.000000f) },
{ (float2)(-0.382683f, -0.923880f) },
{ (float2)(-0.923880f, -0.382683f) },
{ (float2)(-0.707107f, -0.707107f) },
{ (float2)(1.000000f, 0.000000f) },
{ (float2)(0.980785f, -0.195090f) },
{ (float2)(0.831470f, -0.555570f) },
{ (float2)(0.923880f, -0.382683f) },
{ (float2)(0.555570f, -0.831470f) },
{ (float2)(0.707107f, -0.707107f) },
{ (float2)(0.382683f, -0.923880f) },
{ (float2)(0.195090f, -0.980785f) },
{ (float2)(-0.195090f, -0.980785f) },
{ (float2)(0.000000f, -1.000000f) },
{ (float2)(-0.382683f, -0.923880f) },
{ (float2)(-0.555570f, -0.831470f) },
{ (float2)(-0.707107f, -0.707107f) },
{ (float2)(-0.831470f, -0.555570f) },
{ (float2)(-0.980785f, -0.195090f) },
{ (float2)(-0.923880f, -0.382683f) }
};

/* ~\~ end */
/* ~\~ language=OpenCL filename=fftsynth/templates/parity.cl */
/* ~\~ begin <<lit/code-generator.md|fftsynth/templates/parity.cl>>[0] */
inline int parity_2(int i)
{
    int x = MODR(i);

    
    i = DIVR(i);
    x += MODR(i);
    i = DIVR(i);
    x += MODR(i);
    i = DIVR(i);
    x += MODR(i);
    i = DIVR(i);
    x += MODR(i);
    i = DIVR(i);
    x += MODR(i);

    return MODR(x);
}
/* ~\~ end */
/* ~\~ begin <<lit/code-generator.md|fftsynth/templates/parity.cl>>[1] */
#ifdef TESTING
__kernel void test_parity_2(__global const int * x, __global int * y)
{
    int i = get_global_id(0);
    
    y[i] = parity_2(x[i]);
}
#endif // TESTING
/* ~\~ end */
/* ~\~ language=OpenCL filename=fftsynth/templates/transpose.cl */
/* ~\~ begin <<lit/code-generator.md|fftsynth/templates/transpose.cl>>[0] */
inline int transpose_2(int j)
{
    int x = 0;

    
    x = MULR(x) + MODR(j);
    j = DIVR(j);
    x = MULR(x) + MODR(j);
    j = DIVR(j);
    x = MULR(x) + MODR(j);
    j = DIVR(j);
    x = MULR(x) + MODR(j);
    j = DIVR(j);
    x = MULR(x) + MODR(j);
    j = DIVR(j);

    return x;
}
/* ~\~ end */
/* ~\~ begin <<lit/code-generator.md|fftsynth/templates/transpose.cl>>[1] */
#ifdef TESTING
__kernel void test_transpose_2(__global const int * x, __global int * y)
{
    int i = get_global_id(0);
    
    y[i] = transpose_2(x[i]);
}
#endif // TESTING
/* ~\~ end */
/* ~\~ language=OpenCL filename=fftsynth/templates/ipow.cl */
/* ~\~ begin <<lit/code-generator.md|fftsynth/templates/ipow.cl>>[0] */

inline int ipow(int b)
{
    return 1 << b;
}

/* ~\~ end */
/* ~\~ language=OpenCL filename=fftsynth/templates/indices.cl */
/* ~\~ begin <<lit/code-generator.md|fftsynth/templates/indices.cl>>[0] */
inline int comp_idx_2(int i, int k)
{
    int rem = i % ipow(k);
    int base = i - rem;
    return MULR(base) + rem;
}


inline int comp_perm_2(int i, int rem)
{
    int p = parity_2(i);
    return MULR(i) + MODR(rem + 2 - p);
}
/* ~\~ end */
/* ~\~ language=OpenCL filename=fftsynth/templates/fma-codelets.cl */
/* ~\~ begin <<lit/fma-codelets.md|fftsynth/templates/fma-codelets.cl>>[0] */
/* ~\~ begin <<lit/fma-codelets.md|fma-radix2>>[0] */

void fft_2(float8 * restrict s0, float8 * restrict s1, float8 * restrict s0_in, float8 * restrict s1_in, float8 * restrict s0_out, float8 * restrict s1_out, bool first_iteration, bool last_iteration, int cycle, int i0, int i1, int iw)
{
    float8 t0, t1, a, b;
    #ifndef TESTING_RADIX
    __constant float2 *w = W[iw];
    #endif // !TESTING_RADIX
    #ifdef TESTING_RADIX
    float2 w[] = {(float2)(1.0, 0.0)};
    #endif // TESTING_RADIX

    
    switch (cycle) {
        case 1: SWAP(int, i0, i1); break;
    }
    if ( first_iteration )
    {
        t0 = s0_in[i0]; t1 = s1_in[i1];
    }
    else
    {
        t0 = s0[i0]; t1 = s1[i1];
    }
    switch (cycle) {
        case 1: SWAP(float8, t0, t1); break;
    }
    

    a.even = -w[0].odd * t1.odd + t0.even + w[0].even * t1.even;
    a.odd = w[0].odd * t1.even + t0.odd + w[0].even * t1.odd;
    b = 2 * t0 - a;

    t0 = a;
    t1 = b;

    
    switch (cycle) {
        case 1: SWAP(float8, t0, t1); break;
    }
    if ( last_iteration )
    {
        s0_out[i0] = t0; s1_out[i1] = t1;
    }
    else
    {
        s0[i0] = t0; s1[i1] = t1;
    }
    
}

/* ~\~ end */








/* ~\~ begin <<lit/fma-codelets.md|fma-codelet-tests>>[0] */
#ifdef TESTING


__kernel void test_radix_2(__global float8 *x, __global float8 *y, int n)
{
    int i = get_global_id(0) * 2;

    // n is the number of radix2 ffts to perform
    if ( i < 2 * n ) {
        float8 s0 = x[i];
        float8 s1 = x[i + 1];

        fft_2(&s0, &s1, 0, 0, 0, 0);

        y[i] = s0; y[i + 1] = s1;
    }
}







#endif // TESTING
/* ~\~ end */
/* ~\~ end */
/* ~\~ language=OpenCL filename=fftsynth/templates/fma-fft.cl */
/* ~\~ begin <<lit/code-generator.md|fftsynth/templates/fma-fft.cl>>[0] */
void fft_32_ps( float8 * restrict s0, float8 * restrict s1, float8 * restrict s0_in, float8 * restrict s1_in, float8 * restrict s0_out, float8 * restrict s1_out)
{
    int wp = 0;

    for ( uint3_t k = 0; k != 5; ++k )
    {
        int j = (k == 0 ? 0 : ipow(k - 1));

        #pragma ivdep
        for ( uint5_t i = 0; i != 16; ++i )
        {
            int a;
            if ( k != 0 )
            {
                a = comp_idx_2(DIVR(i), k-1);
            }
            else
            {
                a = comp_perm_2(DIVR(i), MODR(i));
            }
            
            fft_2( s0, s1, s0_in, s1_in, s0_out, s1_out, k == 0, k == 4, MODR(i),  a + 0 * j, a + 1 * j, wp);
            
            if ( k != 0 )
            {
                ++wp;
            }
        }
    }
}
/* ~\~ end */

inline void fft_32x32_4(float8 out[restrict 32][32], /*float8 out_upper[restrict 16][32], float8 out_lower[restrict 16][32],*/  const float8 in[restrict 32][32], int sign)
{
#pragma loop_coalesce 2
#pragma ii 1
  for (uint2_t dim = 0; dim < 2; dim ++) {
#pragma ii 1
      for (uint6_t n = 0; n < 32; ++n) {
 		  float8 s0[16], s1[16], s0_in[16], s1_in[16], s0_out[16], s1_out[16];

          for (uint6_t m = 0; m < 32; ++m) {  // m is dimension of fft
              uint6_t xi, yi;
              xi = dim == 0 ? n : m;
              yi = dim == 0 ? m : n;

              int i = transpose_2(m);
              int p = parity_2(i);
			  switch (p)
			  {
				  case 0: s0_in[DIVR(i)] = dim == 0 ? in[xi][yi] : out[xi][yi]; break;
				  case 1: s1_in[DIVR(i)] = dim == 0 ? in[xi][yi] : out[xi][yi]; break;
			  }
          }

          fft_32_ps(s0, s1, s0_in, s1_in, s0_out, s1_out);
          
          for (uint6_t i = 0; i < 32; ++i) {
              int p = parity_2(i);
              uint6_t xi, yi;
              xi = dim == 0 ? n : i;
              yi = dim == 0 ? i : n;
              switch (p) {
                  case 0: out[xi][yi] = s0_out[DIVR(i)]; break;
                  case 1: out[xi][yi] = s1_out[DIVR(i)]; break;
              }
          }
      }
  }
}

