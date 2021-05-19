/* FFT
 * command: python -m fftsynth.generator --radix 4 --depth 3 --fpga --fma
 */
/* ~\~ language=OpenCL filename=fftsynth/templates/preprocessor.cl */
/* ~\~ begin <<lit/code-generator.md|fftsynth/templates/preprocessor.cl>>[0] */
#pragma OPENCL EXTENSION cl_intel_channels : enable
#include <ihc_apint.h>
#ifdef TESTING
channel float2 in_channel, out_channel;
#endif // TESTING
#define SWAP(type, x, y) do { type temp = x; x = y, y = temp; } while ( false );

#define DIVR(x) ((x) >> 2)
#define MODR(x) ((x) & 3)
#define MULR(x) ((x) << 2)

/* ~\~ end */
/* ~\~ language=OpenCL filename=fftsynth/templates/fma-twiddles.cl */
/* ~\~ begin <<lit/fma-codelets.md|fftsynth/templates/fma-twiddles.cl>>[0] */
__constant float2 W[32][2] = {

{ (float2)(1.000000f, 0.000000f), (float2)(1.000000f, 0.000000f) },
{ (float2)(0.923880f, -0.382683f), (float2)(0.707107f, -0.707107f) },
{ (float2)(0.707107f, -0.707107f), (float2)(0.000000f, -1.000000f) },
{ (float2)(0.382683f, -0.923880f), (float2)(-0.707107f, -0.707107f) },
{ (float2)(0.382683f, -0.923880f), (float2)(-0.707107f, -0.707107f) },
{ (float2)(1.000000f, 0.000000f), (float2)(1.000000f, 0.000000f) },
{ (float2)(0.923880f, -0.382683f), (float2)(0.707107f, -0.707107f) },
{ (float2)(0.707107f, -0.707107f), (float2)(0.000000f, -1.000000f) },
{ (float2)(0.707107f, -0.707107f), (float2)(0.000000f, -1.000000f) },
{ (float2)(0.382683f, -0.923880f), (float2)(-0.707107f, -0.707107f) },
{ (float2)(1.000000f, 0.000000f), (float2)(1.000000f, 0.000000f) },
{ (float2)(0.923880f, -0.382683f), (float2)(0.707107f, -0.707107f) },
{ (float2)(0.923880f, -0.382683f), (float2)(0.707107f, -0.707107f) },
{ (float2)(0.707107f, -0.707107f), (float2)(0.000000f, -1.000000f) },
{ (float2)(0.382683f, -0.923880f), (float2)(-0.707107f, -0.707107f) },
{ (float2)(1.000000f, 0.000000f), (float2)(1.000000f, 0.000000f) },
{ (float2)(1.000000f, 0.000000f), (float2)(1.000000f, 0.000000f) },
{ (float2)(0.995185f, -0.098017f), (float2)(0.980785f, -0.195090f) },
{ (float2)(0.980785f, -0.195090f), (float2)(0.923880f, -0.382683f) },
{ (float2)(0.956940f, -0.290285f), (float2)(0.831470f, -0.555570f) },
{ (float2)(0.773010f, -0.634393f), (float2)(0.195090f, -0.980785f) },
{ (float2)(0.923880f, -0.382683f), (float2)(0.707107f, -0.707107f) },
{ (float2)(0.881921f, -0.471397f), (float2)(0.555570f, -0.831470f) },
{ (float2)(0.831470f, -0.555570f), (float2)(0.382683f, -0.923880f) },
{ (float2)(0.555570f, -0.831470f), (float2)(-0.382683f, -0.923880f) },
{ (float2)(0.471397f, -0.881921f), (float2)(-0.555570f, -0.831470f) },
{ (float2)(0.707107f, -0.707107f), (float2)(0.000000f, -1.000000f) },
{ (float2)(0.634393f, -0.773010f), (float2)(-0.195090f, -0.980785f) },
{ (float2)(0.290285f, -0.956940f), (float2)(-0.831470f, -0.555570f) },
{ (float2)(0.195090f, -0.980785f), (float2)(-0.923880f, -0.382683f) },
{ (float2)(0.098017f, -0.995185f), (float2)(-0.980785f, -0.195090f) },
{ (float2)(0.382683f, -0.923880f), (float2)(-0.707107f, -0.707107f) }
};

/* ~\~ end */
/* ~\~ language=OpenCL filename=fftsynth/templates/parity.cl */
/* ~\~ begin <<lit/code-generator.md|fftsynth/templates/parity.cl>>[0] */
inline int parity_4(int i)
{
    int x = MODR(i);

    
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
__kernel void test_parity_4(__global const int * restrict x, __global int * restrict y)
{
    int i = get_global_id(0);

    y[i] = parity_4(x[i]);
}
#endif // TESTING
/* ~\~ end */
/* ~\~ language=OpenCL filename=fftsynth/templates/transpose.cl */
/* ~\~ begin <<lit/code-generator.md|fftsynth/templates/transpose.cl>>[0] */
inline int transpose_4(int j)
{
    int x = 0;

    
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
__kernel void test_transpose_4(__global const int * restrict x, __global int * restrict y)
{
    int i = get_global_id(0);

    y[i] = transpose_4(x[i]);
}
#endif // TESTING
/* ~\~ end */
/* ~\~ language=OpenCL filename=fftsynth/templates/ipow.cl */
/* ~\~ begin <<lit/code-generator.md|fftsynth/templates/ipow.cl>>[0] */

inline int ipow(int b)
{
    return 1 << (2*b);
}

/* ~\~ end */
/* ~\~ language=OpenCL filename=fftsynth/templates/indices.cl */
/* ~\~ begin <<lit/code-generator.md|fftsynth/templates/indices.cl>>[0] */
inline int comp_idx_4(int i, int k)
{
    int rem = i % ipow(k);
    int base = i - rem;
    return MULR(base) + rem;
}


inline int comp_perm_4(int i, int rem)
{
    int p = parity_4(i);
    return MULR(i) + MODR(rem + 4 - p);
}
/* ~\~ end */
/* ~\~ language=OpenCL filename=fftsynth/templates/fma-codelets.cl */
/* ~\~ begin <<lit/fma-codelets.md|fftsynth/templates/fma-codelets.cl>>[0] */
/* ~\~ begin <<lit/fma-codelets.md|fma-radix2>>[0] */

/* ~\~ end */




void fft_4(float2 * restrict s0, float2 * restrict s1, float2 * restrict s2, float2 * restrict s3, float2 * restrict s0_in, float2 * restrict s1_in, float2 * restrict s2_in, float2 * restrict s3_in, float2 * restrict s0_out, float2 * restrict s1_out, float2 * restrict s2_out, float2 * restrict s3_out, bool first_iteration, bool last_iteration, int cycle, int i0, int i1, int i2, int i3, int iw)
{
     float2 t0, t1, t2, t3, a, b, c, d;
    #ifndef TESTING_RADIX
    __constant float2 *w = W[iw];
    #endif // !TESTING_RADIX
    #ifdef TESTING_RADIX
    float2 w[] = {(float2)(1.0, 0.0), (float2)(1.0, 0.0)};
    #endif // TESTING_RADIX

    
    switch (cycle) {
        case 1: SWAP(int, i0, i1); SWAP(int, i2, i3); SWAP(int, i0, i2); break;
        case 2: SWAP(int, i0, i2); SWAP(int, i1, i3); break;
        case 3: SWAP(int, i0, i1); SWAP(int, i1, i3); SWAP(int, i1, i2); break;
    }
    if ( first_iteration )
    {
        t0 = s0_in[i0]; t1 = s1_in[i1]; t2 = s2_in[i2]; t3 = s3_in[i3];
    }
    else
    {
        t0 = s0[i0]; t1 = s1[i1]; t2 = s2[i2]; t3 = s3[i3];
    }
    switch (cycle) {
        case 1: SWAP(float2, t0, t1); SWAP(float2, t1, t2); SWAP(float2, t2, t3); break;
        case 2: SWAP(float2, t0, t2); SWAP(float2, t1, t3); break;
        case 3: SWAP(float2, t2, t3); SWAP(float2, t0, t1); SWAP(float2, t0, t2); break;
    }
    

   // adapted from pedram2013transforming, however
   // some versions of the pedram2013transforming paper, including the one hosted by IEEE and the one hosted here:
   // https://www.cs.utexas.edu/users/flame/pubs/LAC_fft.pdf
   // contains serious errors in the FMA-optimized radix-4 pseudocode
   // finally corrected based on the pseudocode reported in karner1998top

    b.even = t0.even - w[1].even * t2.even + w[1].odd * t2.odd;
    b.odd = t0.odd - w[1].even * t2.odd - w[1].odd * t2.even;
    a.even = 2 * t0.even - b.even;
    a.odd = 2 * t0.odd - b.odd;
    d.even = t1.even - w[1].even * t3.even + w[1].odd * t3.odd;
    d.odd = t1.odd - w[1].even * t3.odd - w[1].odd * t3.even;
    c.even = 2 * t1.even - d.even;
    c.odd = 2 * t1.odd - d.odd;
    t2 = c;
    c.even = a.even - w[0].even * t2.even + w[0].odd * t2.odd;
    c.odd = a.odd - w[0].even * t2.odd - w[0].odd * t2.even;
    t2 = c;
    t0.even = 2 * a.even - c.even;
    t0.odd = 2 * a.odd - c.odd;
    t1 = d;
    d.even = b.even + w[0].even * t1.odd + w[0].odd * t1.even;
    d.odd = b.odd - w[0].even * t1.even + w[0].odd * t1.odd;
    t1 = d;
    t3.even = 2 * b.even - d.even;
    t3.odd = 2 * b.odd - d.odd;

    
    switch (cycle) {
        case 1: SWAP(float2, t2, t3); SWAP(float2, t1, t2); SWAP(float2, t0, t1); break;
        case 2: SWAP(float2, t1, t3); SWAP(float2, t0, t2); break;
        case 3: SWAP(float2, t0, t2); SWAP(float2, t0, t1); SWAP(float2, t2, t3); break;
    }
    if ( last_iteration )
    {
        s0_out[i0] = t0; s1_out[i1] = t1; s2_out[i2] = t2; s3_out[i3] = t3;
    }
    else
    {
        s0[i0] = t0; s1[i1] = t1; s2[i2] = t2; s3[i3] = t3;
    }
    
}





/* ~\~ begin <<lit/fma-codelets.md|fma-codelet-tests>>[0] */
#ifdef TESTING_RADIX






__kernel void test_radix_4(__global float2 *x, __global float2 *y, int n)
{
    int i = get_global_id(0) * 4;

    // n is the number of radix4 ffts to perform
    if (i < 4 * n) {
        float2 s0 = x[i];
        float2 s1 = x[i + 1];
        float2 s2 = x[i + 2];
        float2 s3 = x[i + 3];
        fft_4(&s0, &s1, &s2, &s3, 0, 0, 0, 0, 0, 0);

        y[i] = s0;    y[i + 1] = s1;    y[i + 2] = s2;    y[i + 3] = s3;
    }
}



#endif // TESTING_RADIX
/* ~\~ end */
/* ~\~ end */
/* ~\~ language=OpenCL filename=fftsynth/templates/fma-fft.cl */
/* ~\~ begin <<lit/code-generator.md|fftsynth/templates/fma-fft.cl>>[0] */
void fft_64_ps( float2 * restrict s0, float2 * restrict s1, float2 * restrict s2, float2 * restrict s3, float2 * restrict s0_in, float2 * restrict s1_in, float2 * restrict s2_in, float2 * restrict s3_in, float2 * restrict s0_out, float2 * restrict s1_out, float2 * restrict s2_out, float2 * restrict s3_out)
{
    int wp = 0;

    for ( uint2_t k = 0; k != 3; ++k )
    {
        int j = (k == 0 ? 0 : ipow(k - 1));

        #pragma ivdep
        for ( uint5_t i = 0; i != 16; ++i )
        {
            int a;
            if ( k != 0 )
            {
                a = comp_idx_4(DIVR(i), k-1);
            }
            else
            {
                a = comp_perm_4(DIVR(i), MODR(i));
            }
            
            fft_4( s0, s1, s2, s3, s0_in, s1_in, s2_in, s3_in, s0_out, s1_out, s2_out, s3_out, k == 0, k == 2, MODR(i),  a + 0 * j, a + 1 * j, a + 2 * j, a + 3 * j, wp);
            
            if ( k != 0 )
            {
                ++wp;
            }
        }
    }
}
/* ~\~ end */
/* ~\~ begin <<lit/code-generator.md|fftsynth/templates/fma-fft.cl>>[1] */
#ifdef TESTING
__kernel __attribute__((autorun)) __attribute__((max_global_work_dim(0)))
void fft_64()
{
    while ( true )
    {
    
    float2 s0[16];
    float2 s0_in[16], s0_out[16];
    float2 s1[16];
    float2 s1_in[16], s1_out[16];
    float2 s2[16];
    float2 s2_in[16], s2_out[16];
    float2 s3[16];
    float2 s3_in[16], s3_out[16];

    for ( uint7_t j = 0; j != 64; ++j )
    {
        int i = transpose_4(j);
        int p = parity_4(i);

        float2 x = read_channel_intel(in_channel);
        switch ( p )
        {
            case 0: s0_in[DIVR(i)] = x; break;
            case 1: s1_in[DIVR(i)] = x; break;
            case 2: s2_in[DIVR(i)] = x; break;
            case 3: s3_in[DIVR(i)] = x; break;
        }
    }

    
    fft_64_ps( s0, s1, s2, s3, s0_in, s1_in, s2_in, s3_in, s0_out, s1_out, s2_out, s3_out);
    

    for ( uint7_t i = 0; i != 64; ++i )
    {
        int p = parity_4(i);
        float2 y;
        switch ( p )
        {
            case 0: y = s0_out[DIVR(i)]; break;
            case 1: y = s1_out[DIVR(i)]; break;
            case 2: y = s2_out[DIVR(i)]; break;
            case 3: y = s3_out[DIVR(i)]; break;
        }
        write_channel_intel(out_channel, y);
    }
    }
}
#endif // TESTING
/* ~\~ end */

/* ~\~ language=OpenCL filename=fftsynth/templates/fpga.cl */
/* ~\~ begin <<lit/code-generator.md|fftsynth/templates/fpga.cl>>[0] */
#ifdef TESTING
__kernel __attribute__((max_global_work_dim(0)))
void source(__global const volatile float2 * in, unsigned count)
{
    #pragma ii 1
    for ( unsigned i = 0; i < count; i++ )
    {
        write_channel_intel(in_channel, in[i]);
    }
}

__kernel __attribute__((max_global_work_dim(0)))
void sink(__global float2 *out, unsigned count)
{
    #pragma ii 1
    for ( unsigned i = 0; i < count; i++ )
    {
        out[i] = read_channel_intel(out_channel);
    }
}
#endif // TESTING
/* ~\~ end */

