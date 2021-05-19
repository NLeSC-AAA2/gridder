#include "fft64.cl"

#define NR_RECEIVERS		576
#define	NR_SAMPLES_PER_CHANNEL	3072
#define NR_CHANNELS		64
#define NR_TAPS			16
#define SUBBAND_BANDWIDTH	195312.5f

typedef float2 fcomplex;


inline float2 cmul(float2 a, float2 b)
{
  return (float2)(a.x * b.x - a.y * b.y, a.x * b.y + a.y * b.x);
}


void fft(const float2 *restrict in, float2 *restrict out)
{
  float2 s0[16];
  float2 s0_in[16], s0_out[16];
  float2 s1[16];
  float2 s1_in[16], s1_out[16];
  float2 s2[16];
  float2 s2_in[16], s2_out[16];
  float2 s3[16];
  float2 s3_in[16], s3_out[16];
  
  /* Loop over the number of FFTs to load input */
  for ( uint7_t j = 0; j != 64; ++j )
  {
    int i = transpose_4(j);
    int p = parity_4(i);

    switch ( p )
    {
      case 0: s0_in[DIVR(i)] = in[j]; break;
      case 1: s1_in[DIVR(i)] = in[j]; break;
      case 2: s2_in[DIVR(i)] = in[j]; break;
      case 3: s3_in[DIVR(i)] = in[j]; break;
    }
  }
  /* FFT call */
  fft_64_ps( s0, s1, s2, s3, s0_in, s1_in, s2_in, s3_in, s0_out, s1_out, s2_out, s3_out);
  /* Loop over the number of FFTs to store output */
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
    out[i] = y;
  }
}


typedef char2 InputType[NR_RECEIVERS][NR_SAMPLES_PER_CHANNEL + NR_TAPS - 1][NR_CHANNELS];
typedef float  FIR_FilterWeightsType[NR_CHANNELS][NR_TAPS];
typedef float  BandPassWeightsType[NR_CHANNELS];
typedef struct { float atBegin, afterEnd; } DelaysType[NR_RECEIVERS];
typedef fcomplex OutputType[NR_RECEIVERS][NR_SAMPLES_PER_CHANNEL][NR_CHANNELS]; 


__kernel __attribute__((max_global_work_dim(0)))
void fused
(
  const __global InputType *restrict input,
  __constant FIR_FilterWeightsType *restrict firFilterWeights,
  float subbandFrequency,
  __global const DelaysType *restrict delays,
  __constant BandPassWeightsType *restrict bandPassWeights,
  __global OutputType *restrict output
)
{
  for (unsigned recv = 0; recv < NR_RECEIVERS; recv ++) {
    float2 firFilterOutput[NR_CHANNELS], fftOutput[NR_CHANNELS];
    float2 v[NR_CHANNELS], dv[NR_CHANNELS];
    float2 history[NR_CHANNELS][NR_TAPS];

    float phiBegin = -2 * 3.1415926535f * (*delays)[recv].atBegin;
    float phiEnd = -2 * 3.1415926535f * (*delays)[recv].afterEnd;
    float deltaPhi = (phiEnd - phiBegin) / NR_SAMPLES_PER_CHANNEL;

    for (unsigned ch = 0; ch < NR_CHANNELS; ch ++) {
      float frequency = subbandFrequency - .5f * SUBBAND_BANDWIDTH + ch * (SUBBAND_BANDWIDTH / NR_CHANNELS);
      float myPhiBegin = phiBegin * frequency;
      float myPhiDelta = deltaPhi * frequency;
      v[ch]  = (float2) (native_cos(myPhiBegin), native_sin(myPhiBegin)) * (*bandPassWeights)[ch];
      dv[ch] = (float2) (native_cos(myPhiDelta), native_sin(myPhiDelta));
    }

    for (int time = 1 - NR_TAPS; time < NR_SAMPLES_PER_CHANNEL; time ++) {
#pragma unroll 8
      for (unsigned ch = 0; ch < NR_CHANNELS; ch ++) {
#pragma unroll
	for (unsigned tap = NR_TAPS - 1; tap > 0; tap --)
	  history[ch][tap] = history[ch][tap - 1];

	history[ch][0] = convert_float2((*input)[recv][time + NR_TAPS - 1][ch]);

	float2 sum = 0;

#pragma unroll
	for (unsigned tap = 0; tap < NR_TAPS - 1; tap ++)
	  sum += (*firFilterWeights)[ch][tap] * history[ch][tap];

	firFilterOutput[ch] = sum;
      }

      fft(firFilterOutput, fftOutput);

#pragma unroll 8
      for (unsigned ch = 0; ch < NR_CHANNELS; ch ++) {
	if (time >= 0) {
	  (*output)[recv][time][ch] = cmul(v[ch], fftOutput[ch]);
	  v[ch] = cmul(v[ch], dv[ch]);
	}
      }
    }
  }
}
