#define NR_RECEIVERS		576
#define	NR_SAMPLES_PER_CHANNEL	3072
#define NR_CHANNELS		64
#define NR_TAPS			16
#define SUBBAND_BANDWIDTH	195312.5f


typedef float2 fcomplex;

#define DIR (-1)

typedef float4  float2x2;
typedef float16 float8x2;


inline float2 cmul(float2 a, float2 b)
{
  return (float2)(a.x * b.x - a.y * b.y, a.x * b.y + a.y * b.x);
}


inline float8x2 bitreverse8(float8x2 a)
{
  return (float8x2)(a.s018945CD23AB67EF);
}


void radix2_fwd(float2 *restrict out0, float2 *restrict out1, float2 in0, float2 in1)
{
  *out0 = in0 + in1;
  *out1 = in0 - in1;
}


void radix_8x8_fwd(float2 out[8][8], const float2 in[8][8])
{
#pragma unroll 1
  for (int i = 0; i < 8; i ++) {
    float2 b[8], c[8];

    radix2_fwd(&b[0], &b[4], in[i][0], in[i][4]);
    radix2_fwd(&b[1], &b[5], in[i][1], in[i][5]);
    radix2_fwd(&b[2], &b[6], in[i][2], in[i][6]);
    radix2_fwd(&b[3], &b[7], in[i][3], in[i][7]);

    b[5] = cmul(b[5], (float2)(.70710678118654752440f,(DIR)*.70710678118654752440f));
    b[6] = (float2)(DIR*-b[6].y, DIR*b[6].x);
    b[7] = cmul(b[7], (float2)(-.70710678118654752440f,(DIR)*.70710678118654752440f));

    radix2_fwd(&c[0], &c[2], b[0], b[2]);
    radix2_fwd(&c[1], &c[3], b[1], b[3]);
    radix2_fwd(&c[4], &c[6], b[4], b[6]);
    radix2_fwd(&c[5], &c[7], b[5], b[7]);

    c[3] = (float2)(DIR*-c[3].y, DIR*c[3].x);
    c[7] = (float2)(DIR*-c[7].y, DIR*c[7].x);

    radix2_fwd(&out[0][i], &out[4][i], c[0], c[1]);
    radix2_fwd(&out[2][i], &out[6][i], c[2], c[3]);
    radix2_fwd(&out[1][i], &out[5][i], c[4], c[5]);
    radix2_fwd(&out[3][i], &out[7][i], c[6], c[7]);
  }
}


__constant float2 weights[8][8] = 
{
  {
    (float2)(1.f,0.f),
    (float2)(1.f,0.f),
    (float2)(1.f,0.f),
    (float2)(1.f,0.f),
    (float2)(1.f,0.f),
    (float2)(1.f,0.f),
    (float2)(1.f,0.f),
    (float2)(1.f,0.f),
  },{
    (float2)(1.f,0.f),
    (float2)(0.995184726672f,-0.0980171403296f),
    (float2)(0.980785280403f,-0.195090322016f),
    (float2)(0.956940335732f,-0.290284677254f),
    (float2)(0.923879532511f,-0.382683432365f),
    (float2)(0.881921264348f,-0.471396736826f),
    (float2)(0.831469612303f,-0.55557023302f),
    (float2)(0.773010453363f,-0.634393284164f),
  },{
    (float2)(1.f,0.f),
    (float2)(0.980785280403f,-0.195090322016f),
    (float2)(0.923879532511f,-0.382683432365f),
    (float2)(0.831469612303f,-0.55557023302f),
    (float2)(0.707106781187f,-0.707106781187f),
    (float2)(0.55557023302f,-0.831469612303f),
    (float2)(0.382683432365f,-0.923879532511f),
    (float2)(0.195090322016f,-0.980785280403f),
  },{
    (float2)(1.f,0.f),
    (float2)(0.956940335732f,-0.290284677254f),
    (float2)(0.831469612303f,-0.55557023302f),
    (float2)(0.634393284164f,-0.773010453363f),
    (float2)(0.382683432365f,-0.923879532511f),
    (float2)(0.0980171403296f,-0.995184726672f),
    (float2)(-0.195090322016f,-0.980785280403f),
    (float2)(-0.471396736826f,-0.881921264348f),
  },{
    (float2)(1.f,0.f),
    (float2)(0.923879532511f,-0.382683432365f),
    (float2)(0.707106781187f,-0.707106781187f),
    (float2)(0.382683432365f,-0.923879532511f),
    (float2)(0.f,-1.f),
    (float2)(-0.382683432365f,-0.923879532511f),
    (float2)(-0.707106781187f,-0.707106781187f),
    (float2)(-0.923879532511f,-0.382683432365f),
  },{
    (float2)(1.f,0.f),
    (float2)(0.881921264348f,-0.471396736826f),
    (float2)(0.55557023302f,-0.831469612303f),
    (float2)(0.0980171403296f,-0.995184726672f),
    (float2)(-0.382683432365f,-0.923879532511f),
    (float2)(-0.773010453363f,-0.634393284164f),
    (float2)(-0.980785280403f,-0.195090322016f),
    (float2)(-0.956940335732f,0.290284677254f),
  },{
    (float2)(1.f,0.f),
    (float2)(0.831469612303f,-0.55557023302f),
    (float2)(0.382683432365f,-0.923879532511f),
    (float2)(-0.195090322016f,-0.980785280403f),
    (float2)(-0.707106781187f,-0.707106781187f),
    (float2)(-0.980785280403f,-0.195090322016f),
    (float2)(-0.923879532511f,0.382683432365f),
    (float2)(-0.55557023302f,0.831469612303f),
  },{
    (float2)(1.f,0.f),
    (float2)(0.773010453363f,-0.634393284164f),
    (float2)(0.195090322016f,-0.980785280403f),
    (float2)(-0.471396736826f,-0.881921264348f),
    (float2)(-0.923879532511f,-0.382683432365f),
    (float2)(-0.956940335732f,0.290284677254f),
    (float2)(-0.55557023302f,0.831469612303f),
    (float2)(0.0980171403296f,0.995184726672f),
  }
};


void fft(const float2 *restrict in, float2 *restrict out)
{
  float2 a[8][8], b[8][8], c[8][8], d[8][8];

#pragma unroll
  for (int i = 0; i < 8; i ++)
#pragma unroll
    for (int j = 0; j < 8; j ++)
      a[j][i] = in[8*i+j];

  radix_8x8_fwd(b, a);

#pragma unroll
  for (int i = 0; i < 8; i ++)
#pragma unroll
    for (int j = 0; j < 8; j ++)
      c[i][j] = cmul(b[i][j], weights[i][j]);

  radix_8x8_fwd(d, c);

#pragma unroll
  for (int i = 0; i < 8; i ++)
#pragma unroll
    for (int j = 0; j < 8; j ++)
      out[8 * i + j] = d[i][j];
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
