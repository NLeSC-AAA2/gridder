#include <ihc_apint.h>

#define SUBGRID_SIZE			32
#define NR_TIMESTEPS			128
#define NR_CHANNELS			16
#define NR_STATIONS			64
#define NR_POLARIZATIONS		4
#define UNROLL_PIXELS_FACTOR		4
#define UNROLL_VISIBILITIES_FACTOR	4

#define NR_PIXELS			((SUBGRID_SIZE) * (SUBGRID_SIZE))

#define REAL				0
#define IMAG				1
#define COMPLEX				2

#define MAX(A,B)			((A)>(B)?(A):(B))
#define NEXT_POWER_OF_TWO(N)		((((N)-1)|((N-1)>>1)|((N)-1)>>2|((N)-1)>>3|((N)-1)>>4|((N)-1)>>5|((N)-1)>>6|((N)-1)>>7|((N)-1)>>8)+1)



inline float2 cmul(float2 a, float2 b)
{
  return (float2) (a.x * b.x + -a.y * b.y, a.x * b.y + a.y * b.x);
}


inline float2 conj(float2 a)
{
  return (float2) (a.x, -a.y);
}


union matrix
{
  float8 f8;
  float2 f2[2][2];
};


inline union matrix matmul(union matrix a, union matrix b)
{
  //float2 a[2][2] = {{ matrixA.s01, matrixA.s23 }, { matrixA.s45, matrixA.s67 }};
  //float2 b[2][2] = {{ matrixB.s01, matrixB.s23 }, { matrixB.s45, matrixB.s67 }};
  //float2 c[2][2];
  union matrix c;

#pragma loop_coalesce
#pragma ii 1
  for (uint2_t y = 0; y < 2; y ++) {
#pragma ii 1
    for (uint2_t x = 0; x < 2; x ++) {
      float2 sum = 0;

#pragma unroll
      for (uint2_t k = 0; k < 2; k ++) {
	sum += cmul(a.f2[y][k], b.f2[k][x]);
      }

      c.f2[y][x] = sum;
    }
  }

  //return (float8) (c[0][0], c[0][1], c[1][0], c[1][1]);
  return c;
}


inline union matrix conjugate(union matrix matrix)
{
  matrix.f8.odd = -matrix.f8.odd;
  return matrix;
}


inline union matrix transpose(union matrix matrix)
{
  matrix.f8 = matrix.f8.s01452367;
  return matrix;
}


inline union matrix hermitian(union matrix matrix)
{
  return transpose(conjugate(matrix));
}


#if defined USE_SIN_COS_LOOKUP_TABLE

#define NR_BITS     11
#define NR_ENTRIES  (1 << (NR_BITS))
#define BIT(A,N)    (((A) >> (N)) & 1)

__constant uint60_t cosisin_table[NR_ENTRIES] = {
  0xfe0000000000000,
  0xfdffffff9c90fdb,
  0xfdffffefa490fda,
  0xfdffffd7a96cbe2,
  0xfdffffb3ac90fd6,
  0xfdffff87afb53c7,
  0xfdffff53b16cbdb,
  0xfdffff13b2fedd1,
  0xfdfffec7b490fc6,
  0xfdfffe73b6231b9,
  0xfdfffe17b7b53a9,
  0xfdfffdafb8a3acb,
  0xfdfffd3bb96cbc1,
  0xfdfffcc3ba35cb6,
  0xfdfffc3bbafeda8,
  0xfdfffbabbbc7e99,
  0xfdfffb13bc90f88,
  0xfdfffa6fbd5a075,
  0xfdfff9c3be23160,
  0xfdfff90fbeec24a,
  0xfdfff84fbfb5330,
  0xfdfff783c03f20a,
  0xfdfff6afc0a3a7b,
  0xfdfff5cfc1082ea,
  0xfdfff4e7c16cb58,
  0xfdfff3f7c1d13c5,
  0xfdfff2fbc235c31,
  0xfdfff1f7c29a499,
  0xfdfff0e7c2fed01,
  0xfdffefcbc363568,
  0xfdffeeabc3c7dcc,
  0xfdffed7bc42c62f,
  0xfdffec47c490e90,
  0xfdffeb07c4f56ee,
  0xfdffe9bbc559f4b,
  0xfdffe867c5be7a6,
  0xfdffe707c622fff,
  0xfdffe59fc687856,
  0xfdffe42fc6ec0aa,
  0xfdffe2b3c7508fb,
  0xfdffe12bc7b514b,
  0xfdffdf9bc80cccc,
  0xfdffde03c83f0f1,
  0xfdffdc5fc871516,
  0xfdffdab3c8a3938,
  0xfdffd8fbc8d5d5a,
  0xfdffd73bc90817a,
  0xfdffd56fc93a599,
  0xfdffd39bc96c9b6,
  0xfdffd1bbc99edd1,
  0xfdffcfd3c9d11ec,
  0xfdffcddfca03605,
  0xfdffcbe3ca35a1c,
  0xfdffc9dfca67e32,
  0xfdffc7cfca9a245,
  0xfdffc5b3cacc658,
  0xfdffc38fcafea69,
  0xfdffc163cb30e78,
  0xfdffbf2bcb63285,
  0xfdffbcebcb95692,
  0xfdffba9fcbc7a9b,
  0xfdffb84bcbf9ea3,
  0xfdffb5ebcc2c2a9,
  0xfdffb383cc5e6ad,
  0xfdffb10fcc90ab0,
  0xfdffae93ccc2eb0,
  0xfdffac0bccf52ae,
  0xfdffa97bcd276ab,
  0xfdffa6e3cd59aa5,
  0xfdffa43fcd8be9e,
  0xfdffa18fcdbe294,
  0xfdff9ed7cdf0689,
  0xfdff9c17ce22a7a,
  0xfdff994bce54e6a,
  0xfdff9677ce87258,
  0xfdff9397ceb9643,
  0xfdff90afceeba2d,
  0xfdff8dbbcf1de13,
  0xfdff8abfcf501f7,
  0xfdff87b7cf825da,
  0xfdff84a7cfb49b9,
  0xfdff818bcfe6d97,
  0xfdff7e67d00c8b9,
  0xfdff7b3bd025aa5,
  0xfdff7803d03ec90,
  0xfdff74bfd057e79,
  0xfdff7173d071062,
  0xfdff6e1fd08a249,
  0xfdff6abfd0a342f,
  0xfdff6757d0bc613,
  0xfdff63e3d0d57f7,
  0xfdff6067d0ee9d8,
  0xfdff5cdfd107bb8,
  0xfdff594fd120d97,
  0xfdff55b3d139f75,
  0xfdff520fd153151,
  0xfdff4e5fd16c32c,
  0xfdff4aa7d185505,
  0xfdff46e7d19e6dc,
  0xfdff431bd1b78b3,
  0xfdff3f47d1d0a88,
  0xfdff3b67d1e9c5b,
  0xfdff377bd202e2c,
  0xfdff3387d21bffd,
  0xfdff2f8bd2351cc,
  0xfdff2b83d24e398,
  0xfdff2773d267564,
  0xfdff2357d28072e,
  0xfdff1f33d2998f6,
  0xfdff1b07d2b2abd,
  0xfdff16cbd2cbc82,
  0xfdff128bd2e4e45,
  0xfdff0e3fd2fe006,
  0xfdff09e7d3171c6,
  0xfdff0587d330385,
  0xfdff011fd349541,
  0xfdfefcabd3626fc,
  0xfdfef82fd37b8b5,
  0xfdfef3a7d394a6d,
  0xfdfeef17d3adc21,
  0xfdfeea7bd3c6dd5,
  0xfdfee5d7d3dff87,
  0xfdfee127d3f9137,
  0xfdfedc6fd4122e5,
  0xfdfed7abd42b492,
  0xfdfed2dfd44463c,
  0xfdfece0bd45d7e4,
  0xfdfec92bd47698b,
  0xfdfec43fd48fb30,
  0xfdfebf4bd4a8cd2,
  0xfdfeba4fd4c1e73,
  0xfdfeb547d4db012,
  0xfdfeb037d4f41ae,
  0xfdfeab1bd50d349,
  0xfdfea5f7d5264e2,
  0xfdfea0c7d53f67a,
  0xfdfe9b8fd55880e,
  0xfdfe964bd5719a1,
  0xfdfe90ffd58ab32,
  0xfdfe8ba7d5a3cbf,
  0xfdfe8647d5bce4c,
  0xfdfe80dfd5d5fd7,
  0xfdfe7b6bd5ef15f,
  0xfdfe75ebd6082e4,
  0xfdfe7063d621469,
  0xfdfe6ad3d63a5eb,
  0xfdfe6537d653769,
  0xfdfe5f93d66c8e7,
  0xfdfe59e3d685a62,
  0xfdfe542bd69ebdb,
  0xfdfe4e67d6b7d51,
  0xfdfe489bd6d0ec6,
  0xfdfe42c7d6ea038,
  0xfdfe3ce7d7031a7,
  0xfdfe36fbd71c315,
  0xfdfe3107d735480,
  0xfdfe2b0bd74e5e8,
  0xfdfe2503d76774f,
  0xfdfe1eefd7808b3,
  0xfdfe18d3d799a15,
  0xfdfe12afd7b2b74,
  0xfdfe0c7fd7cbcd1,
  0xfdfe0647d7e4e2b,
  0xfdfe0003d7fdf82,
  0xfdfdf9b7d80b86c,
  0xfdfdf363d818116,
  0xfdfded03d8249be,
  0xfdfde697d831265,
  0xfdfde023d83db0a,
  0xfdfdd9a7d84a3af,
  0xfdfdd31fd856c52,
  0xfdfdcc8bd8634f4,
  0xfdfdc5efd86fd95,
  0xfdfdbf4bd87c634,
  0xfdfdb89bd888ed2,
  0xfdfdb1e3d89576e,
  0xfdfdab1fd8a2009,
  0xfdfda453d8ae8a3,
  0xfdfd9d7fd8bb13c,
  0xfdfd969fd8c79d4,
  0xfdfd8fb3d8d426a,
  0xfdfd88bfd8e0aff,
  0xfdfd81bfd8ed391,
  0xfdfd7abbd8f9c23,
  0xfdfd73a7d9064b4,
  0xfdfd6c8bd912d43,
  0xfdfd6567d91f5d1,
  0xfdfd5e37d92be5d,
  0xfdfd56ffd9386e8,
  0xfdfd4fbbd944f71,
  0xfdfd486fd9517f9,
  0xfdfd4117d95e07f,
  0xfdfd39b7d96a905,
  0xfdfd324fd977188,
  0xfdfd2adbd983a0b,
  0xfdfd235bd99028c,
  0xfdfd1bd3d99cb0a,
  0xfdfd1443d9a9388,
  0xfdfd0ca7d9b5c04,
  0xfdfd0503d9c247f,
  0xfdfcfd53d9cecf9,
  0xfdfcf59bd9db571,
  0xfdfcedd7d9e7de7,
  0xfdfce60bd9f465b,
  0xfdfcde37da00ece,
  0xfdfcd653da0d740,
  0xfdfcce6bda19fb0,
  0xfdfcc677da2681f,
  0xfdfcbe77da3308c,
  0xfdfcb673da3f8f7,
  0xfdfcae5fda4c161,
  0xfdfca643da589c9,
  0xfdfc9e1fda65230,
  0xfdfc95efda71a95,
  0xfdfc8db7da7e2f9,
  0xfdfc8573da8ab5b,
  0xfdfc7d27da973ba,
  0xfdfc74d3daa3c19,
  0xfdfc6c73dab0476,
  0xfdfc6407dabccd1,
  0xfdfc5b93dac952b,
  0xfdfc5317dad5d83,
  0xfdfc4a8fdae25d9,
  0xfdfc41ffdaeee2d,
  0xfdfc3963dafb680,
  0xfdfc30bfdb07ed1,
  0xfdfc280fdb14721,
  0xfdfc1f57db20f6f,
  0xfdfc1693db2d7bb,
  0xfdfc0dc7db3a006,
  0xfdfc04f3db4684d,
  0xfdfbfc13db53095,
  0xfdfbf327db5f8da,
  0xfdfbea33db6c11d,
  0xfdfbe137db7895f,
  0xfdfbd82fdb8519f,
  0xfdfbcf1fdb919dd,
  0xfdfbc603db9e219,
  0xfdfbbcdfdbaaa53,
  0xfdfbb3afdbb728c,
  0xfdfbaa77dbc3ac3,
  0xfdfba133dbd02f9,
  0xfdfb97e7dbdcb2c,
  0xfdfb8e93dbe935e,
  0xfdfb8533dbf5b8c,
  0xfdfb7bcbdc023ba,
  0xfdfb7257dc0ebe6,
  0xfdfb68d7dc1b410,
  0xfdfb5f53dc27c39,
  0xfdfb55bfdc3445f,
  0xfdfb4c27dc40c84,
  0xfdfb4283dc4d4a5,
  0xfdfb38d3dc59cc6,
  0xfdfb2f1bdc664e5,
  0xfdfb2557dc72d02,
  0xfdfb1b8fdc7f51d,
  0xfdfb11b7dc8bd36,
  0xfdfb07d7dc9854c,
  0xfdfafdefdca4d62,
  0xfdfaf3fbdcb1575,
  0xfdfae9ffdcbdd86,
  0xfdfadff7dcca596,
  0xfdfad5e7dcd6da3,
  0xfdfacbcbdce35ae,
  0xfdfac1a7dcefdb7,
  0xfdfab77bdcfc5be,
  0xfdfaad43dd08dc4,
  0xfdfaa303dd155c7,
  0xfdfa98b7dd21dc9,
  0xfdfa8e5fdd2e5c8,
  0xfdfa8403dd3adc5,
  0xfdfa7997dd475c0,
  0xfdfa6f27dd53db9,
  0xfdfa64a7dd605b0,
  0xfdfa5a23dd6cda5,
  0xfdfa4f93dd79598,
  0xfdfa44f7dd85d89,
  0xfdfa3a53dd92578,
  0xfdfa2fa7dd9ed64,
  0xfdfa24efddab54f,
  0xfdfa1a2fddb7d37,
  0xfdfa0f63ddc451e,
  0xfdfa048fddd0d02,
  0xfdf9f9afdddd4e4,
  0xfdf9eec7dde9cc4,
  0xfdf9e3d3ddf64a1,
  0xfdf9d8d7de02c7d,
  0xfdf9cdd3de0f457,
  0xfdf9c2c3de1bc2e,
  0xfdf9b7abde28404,
  0xfdf9ac87de34bd7,
  0xfdf9a157de413a8,
  0xfdf99623de4db75,
  0xfdf98adfde5a342,
  0xfdf97f97de66b0c,
  0xfdf97443de732d4,
  0xfdf968e3de7fa9a,
  0xfdf95d7bde8c25e,
  0xfdf9520bde98a1f,
  0xfdf9468fdea51dd,
  0xfdf93b07deb199a,
  0xfdf92f7bdebe155,
  0xfdf923dfdeca90d,
  0xfdf9183fded70c3,
  0xfdf90c8fdee3877,
  0xfdf900dbdef0027,
  0xfdf8f51bdefc7d6,
  0xfdf8e94fdf08f83,
  0xfdf8dd7bdf1572e,
  0xfdf8d19fdf21ed6,
  0xfdf8c5b7df2e67c,
  0xfdf8b9c7df3ae1f,
  0xfdf8adcbdf475bf,
  0xfdf8a1c7df53d5e,
  0xfdf895b7df604fb,
  0xfdf8899fdf6cc95,
  0xfdf87d7bdf7942d,
  0xfdf8714fdf85bc2,
  0xfdf8651bdf92355,
  0xfdf858dbdf9eae4,
  0xfdf84c93dfab272,
  0xfdf8403fdfb79fe,
  0xfdf833dfdfc4187,
  0xfdf8277bdfd090e,
  0xfdf81b0bdfdd093,
  0xfdf80e8fdfe9815,
  0xfdf8020bdff5f93,
  0xfdf7f57be001388,
  0xfdf7e8e3e007745,
  0xfdf7dc43e00db01,
  0xfdf7cf97e013ebc,
  0xfdf7c2e3e01a276,
  0xfdf7b623e02062e,
  0xfdf7a95be0269e5,
  0xfdf79c87e02cd9b,
  0xfdf78fabe03314f,
  0xfdf782c7e039502,
  0xfdf775d7e03f8b4,
  0xfdf768dbe045c65,
  0xfdf75bd7e04c014,
  0xfdf74ecbe0523c2,
  0xfdf741b3e05876f,
  0xfdf73493e05eb1b,
  0xfdf7276be064ec5,
  0xfdf71a33e06b26e,
  0xfdf70cf7e071616,
  0xfdf6ffafe0779bc,
  0xfdf6f25be07dd61,
  0xfdf6e503e084105,
  0xfdf6d79be08a4a8,
  0xfdf6ca2be090848,
  0xfdf6bcb3e096be8,
  0xfdf6af33e09cf86,
  0xfdf6a1a3e0a3323,
  0xfdf6940fe0a96bf,
  0xfdf6866fe0afa5a,
  0xfdf678c3e0b5df3,
  0xfdf66b13e0bc18b,
  0xfdf65d53e0c2521,
  0xfdf64f8be0c88b6,
  0xfdf641bbe0cec4a,
  0xfdf633e3e0d4fdd,
  0xfdf625fbe0db36e,
  0xfdf6180fe0e16fe,
  0xfdf60a17e0e7a8b,
  0xfdf5fc13e0ede18,
  0xfdf5ee0be0f41a3,
  0xfdf5dff3e0fa52e,
  0xfdf5d1d3e1008b7,
  0xfdf5c3abe106c3e,
  0xfdf5b57be10cfc4,
  0xfdf5a73be113349,
  0xfdf598f7e1196cc,
  0xfdf58aa7e11fa4e,
  0xfdf57c4fe125dce,
  0xfdf56debe12c14d,
  0xfdf55f7be1324cb,
  0xfdf55107e138846,
  0xfdf54283e13ebc1,
  0xfdf533fbe144f3a,
  0xfdf52567e14b2b2,
  0xfdf516c7e151629,
  0xfdf5081fe15799e,
  0xfdf4f96fe15dd11,
  0xfdf4eab3e164083,
  0xfdf4dbefe16a3f4,
  0xfdf4cd1fe170763,
  0xfdf4be47e176ad1,
  0xfdf4af63e17ce3e,
  0xfdf4a077e1831a9,
  0xfdf49183e189512,
  0xfdf48283e18f879,
  0xfdf47377e195bdf,
  0xfdf46467e19bf44,
  0xfdf45547e1a22a8,
  0xfdf44623e1a860a,
  0xfdf436f3e1ae96b,
  0xfdf427b7e1b4cca,
  0xfdf41873e1bb027,
  0xfdf40927e1c1383,
  0xfdf3f9cfe1c76de,
  0xfdf3ea6be1cda37,
  0xfdf3db03e1d3d8e,
  0xfdf3cb8be1da0e4,
  0xfdf3bc0fe1e0439,
  0xfdf3ac87e1e678b,
  0xfdf39cf3e1ecadc,
  0xfdf38d57e1f2e2c,
  0xfdf37db3e1f917b,
  0xfdf36e03e1ff4c7,
  0xfdf35e4be205813,
  0xfdf34e87e20bb5c,
  0xfdf33ebbe211ea5,
  0xfdf32ee3e2181eb,
  0xfdf31f03e21e530,
  0xfdf30f1be224874,
  0xfdf2ff27e22abb6,
  0xfdf2ef2be230ef6,
  0xfdf2df23e237234,
  0xfdf2cf13e23d571,
  0xfdf2bef7e2438ad,
  0xfdf2aed3e249be7,
  0xfdf29ea7e24ff20,
  0xfdf28e6fe256257,
  0xfdf27e2be25c58c,
  0xfdf26de3e2628c0,
  0xfdf25d8be268bf2,
  0xfdf24d2fe26ef22,
  0xfdf23cc7e275251,
  0xfdf22c53e27b57f,
  0xfdf21bd7e2818aa,
  0xfdf20b53e287bd4,
  0xfdf1fac3e28defc,
  0xfdf1ea2be294223,
  0xfdf1d987e29a548,
  0xfdf1c8dbe2a086b,
  0xfdf1b823e2a6b8d,
  0xfdf1a763e2acead,
  0xfdf1969be2b31cc,
  0xfdf185c7e2b94e9,
  0xfdf174ebe2bf804,
  0xfdf16403e2c5b1e,
  0xfdf15313e2cbe36,
  0xfdf14217e2d214c,
  0xfdf13113e2d8461,
  0xfdf12007e2de774,
  0xfdf10eefe2e4a84,
  0xfdf0fdcfe2ead94,
  0xfdf0eca3e2f10a2,
  0xfdf0db6fe2f73ae,
  0xfdf0ca2fe2fd6b9,
  0xfdf0b8e7e3039c2,
  0xfdf0a797e309cc9,
  0xfdf0963be30ffce,
  0xfdf084d3e3162d2,
  0xfdf07367e31c5d4,
  0xfdf061ebe3228d4,
  0xfdf0506be328bd3,
  0xfdf03edfe32eed0,
  0xfdf02d47e3351cb,
  0xfdf01ba7e33b4c4,
  0xfdf009ffe3417bb,
  0xfdeff84be347ab1,
  0xfdefe68fe34dda6,
  0xfdefd4c7e354098,
  0xfdefc2f7e35a389,
  0xfdefb11fe360678,
  0xfdef9f3be366965,
  0xfdef8d4fe36cc51,
  0xfdef7b57e372f3a,
  0xfdef6957e379222,
  0xfdef574be37f508,
  0xfdef4537e3857ed,
  0xfdef3317e38bacf,
  0xfdef20f3e391daf,
  0xfdef0ebfe39808e,
  0xfdeefc83e39e36c,
  0xfdeeea3fe3a4647,
  0xfdeed7f3e3aa921,
  0xfdeec59be3b0bf9,
  0xfdeeb337e3b6ecf,
  0xfdeea0cbe3bd1a3,
  0xfdee8e57e3c3476,
  0xfdee7bd7e3c9746,
  0xfdee694fe3cfa15,
  0xfdee56bbe3d5ce2,
  0xfdee441fe3dbfae,
  0xfdee317be3e2276,
  0xfdee1ecbe3e853d,
  0xfdee0c13e3ee803,
  0xfdedf94fe3f4ac7,
  0xfdede683e3fad89,
  0xfdedd3abe401049,
  0xfdedc0cbe407308,
  0xfdedade3e40d5c4,
  0xfded9aefe41387f,
  0xfded87f3e419b37,
  0xfded74ebe41fdee,
  0xfded61dbe4260a3,
  0xfded4ec3e42c357,
  0xfded3b9fe432608,
  0xfded286fe4388b6,
  0xfded153be43eb64,
  0xfded01f7e444e0f,
  0xfdeceeafe44b0b9,
  0xfdecdb5be451361,
  0xfdecc7fbe457607,
  0xfdecb493e45d8ab,
  0xfdeca123e463b4d,
  0xfdec8da7e469ded,
  0xfdec7a23e47008b,
  0xfdec6693e476328,
  0xfdec52fbe47c5c2,
  0xfdec3f5be48285b,
  0xfdec2bafe488af0,
  0xfdec17fbe48ed85,
  0xfdec043be495018,
  0xfdebf073e49b2a9,
  0xfdebdca3e4a1537,
  0xfdebc8c7e4a77c4,
  0xfdebb4dfe4ada4f,
  0xfdeba0f3e4b3cd8,
  0xfdeb8cf7e4b9f5f,
  0xfdeb78f7e4c01e4,
  0xfdeb64ebe4c6467,
  0xfdeb50d3e4cc6e9,
  0xfdeb3cb3e4d2968,
  0xfdeb288be4d8be5,
  0xfdeb145be4dee5f,
  0xfdeb001be4e50d8,
  0xfdeaebd7e4eb350,
  0xfdead787e4f15c5,
  0xfdeac32fe4f7838,
  0xfdeaaecbe4fdaa9,
  0xfdea9a5fe503d19,
  0xfdea85e7e509f86,
  0xfdea7167e5101f1,
  0xfdea5cdfe51645a,
  0xfdea484be51c6c1,
  0xfdea33abe522927,
  0xfdea1f07e528b8a,
  0xfdea0a57e52edeb,
  0xfde9f59be535049,
  0xfde9e0d7e53b2a6,
  0xfde9cc0be541501,
  0xfde9b733e54775a,
  0xfde9a253e54d9b1,
  0xfde98d67e553c06,
  0xfde97873e559e58,
  0xfde96377e5600a9,
  0xfde94e6fe5662f8,
  0xfde9395fe56c544,
  0xfde92443e57278f,
  0xfde90f1fe5789d7,
  0xfde8f9f3e57ec1e,
  0xfde8e4bbe584e61,
  0xfde8cf7be58b0a3,
  0xfde8ba2fe5912e4,
  0xfde8a4dbe597522,
  0xfde88f7be59d75e,
  0xfde87a13e5a3997,
  0xfde864a3e5a9bcf,
  0xfde84f27e5afe05,
  0xfde839a3e5b6038,
  0xfde82413e5bc26a,
  0xfde80e7be5c2499,
  0xfde7f8dbe5c86c6,
  0xfde7e32fe5ce8f1,
  0xfde7cd7be5d4b1a,
  0xfde7b7bbe5dad40,
  0xfde7a1f3e5e0f65,
  0xfde78c23e5e7187,
  0xfde77647e5ed3a7,
  0xfde76063e5f35c6,
  0xfde74a73e5f97e2,
  0xfde7347be5ff9fc,
  0xfde71e77e605c13,
  0xfde7086be60be29,
  0xfde6f257e61203c,
  0xfde6dc37e61824e,
  0xfde6c60fe61e45d,
  0xfde6afdfe62466a,
  0xfde699a3e62a874,
  0xfde6835fe630a7c,
  0xfde66d0fe636c82,
  0xfde656b7e63ce86,
  0xfde64053e643088,
  0xfde629e7e649288,
  0xfde61373e64f486,
  0xfde5fcf3e655681,
  0xfde5e66be65b87a,
  0xfde5cfd7e661a71,
  0xfde5b93be667c66,
  0xfde5a297e66de58,
  0xfde58be7e674049,
  0xfde5752fe67a237,
  0xfde55e6be680422,
  0xfde5479fe68660b,
  0xfde530cbe68c7f2,
  0xfde519ebe6929d7,
  0xfde50303e698bba,
  0xfde4ec0fe69ed9b,
  0xfde4d513e6a4f79,
  0xfde4be0fe6ab155,
  0xfde4a6ffe6b132f,
  0xfde48fe7e6b7507,
  0xfde478c3e6bd6dc,
  0xfde46197e6c38af,
  0xfde44a63e6c9a80,
  0xfde43323e6cfc4e,
  0xfde41bdbe6d5e19,
  0xfde40487e6dbfe3,
  0xfde3ed2be6e21ab,
  0xfde3d5c7e6e8370,
  0xfde3be57e6ee533,
  0xfde3a6dfe6f46f4,
  0xfde38f5be6fa8b2,
  0xfde377cfe700a6e,
  0xfde36037e706c28,
  0xfde3489be70cddf,
  0xfde330efe712f94,
  0xfde3193fe719147,
  0xfde30183e71f2f8,
  0xfde2e9bbe7254a6,
  0xfde2d1efe72b650,
  0xfde2ba13e7317fa,
  0xfde2a233e7379a1,
  0xfde28a47e73db46,
  0xfde2724fe743ce8,
  0xfde25a53e749e88,
  0xfde24247e750026,
  0xfde22a37e7561c1,
  0xfde2121be75c35a,
  0xfde1f9f7e7624f1,
  0xfde1e1c7e768685,
  0xfde1c98fe76e817,
  0xfde1b14be7749a6,
  0xfde198ffe77ab34,
  0xfde180abe780cbd,
  0xfde1684be786e46,
  0xfde14fe3e78cfcc,
  0xfde1376fe79314f,
  0xfde11ef3e7992d0,
  0xfde1066fe79f44f,
  0xfde0eddfe7a55cb,
  0xfde0d547e7ab745,
  0xfde0bca7e7b18bd,
  0xfde0a3fbe7b7a32,
  0xfde08b47e7bdba4,
  0xfde07287e7c3d15,
  0xfde059bfe7c9e82,
  0xfde040efe7cffed,
  0xfde02813e7d6156,
  0xfde00f2be7dc2bc,
  0xfddff63fe7e2420,
  0xfddfdd47e7e8582,
  0xfddfc443e7ee6e1,
  0xfddfab3be7f483d,
  0xfddf9223e7fa998,
  0xfddf7907e800578,
  0xfddf5fdfe803622,
  0xfddf46afe8066cc,
  0xfddf2d73e809774,
  0xfddf142fe80c81b,
  0xfddefadfe80f8c0,
  0xfddee187e812965,
  0xfddec827e815a08,
  0xfddeaebbe818aaa,
  0xfdde9547e81bb4b,
  0xfdde7bcbe81ebea,
  0xfdde6243e821c88,
  0xfdde48b3e824d25,
  0xfdde2f17e827dc0,
  0xfdde1573e82ae5a,
  0xfdddfbc7e82def4,
  0xfddde20fe830f8b,
  0xfdddc84fe834022,
  0xfdddae83e8370b7,
  0xfddd94afe83a14c,
  0xfddd7ad3e83d1de,
  0xfddd60ebe840270,
  0xfddd46fbe843300,
  0xfddd2d03e84638f,
  0xfddd12ffe84941d,
  0xfddcf8f3e84c4aa,
  0xfddcdedbe84f535,
  0xfddcc4bbe8525bf,
  0xfddcaa93e855647,
  0xfddc905fe8586cf,
  0xfddc7623e85b755,
  0xfddc5bdbe85e7d9,
  0xfddc418be86185d,
  0xfddc2733e8648df,
  0xfddc0cd3e867960,
  0xfddbf263e86a9e0,
  0xfddbd7efe86da5e,
  0xfddbbd6fe870adb,
  0xfddba2e7e873b57,
  0xfddb8857e876bd1,
  0xfddb6dbbe879c49,
  0xfddb5313e87ccc1,
  0xfddb3867e87fd37,
  0xfddb1dafe882dac,
  0xfddb02ebe885e20,
  0xfddae81fe888e93,
  0xfddacd4be88bf04,
  0xfddab26be88ef74,
  0xfdda9783e891fe2,
  0xfdda7c93e895050,
  0xfdda6197e8980bc,
  0xfdda4693e89b126,
  0xfdda2b87e89e190,
  0xfdda106fe8a11f7,
  0xfdd9f54fe8a425e,
  0xfdd9da23e8a72c3,
  0xfdd9beefe8aa327,
  0xfdd9a3b3e8ad38a,
  0xfdd9886be8b03eb,
  0xfdd96d1be8b344b,
  0xfdd951bfe8b64aa,
  0xfdd9365fe8b9507,
  0xfdd91aefe8bc563,
  0xfdd8ff7be8bf5bd,
  0xfdd8e3fbe8c2617,
  0xfdd8c86fe8c566f,
  0xfdd8acdfe8c86c5,
  0xfdd89143e8cb71a,
  0xfdd8759be8ce76d,
  0xfdd859ebe8d17bf,
  0xfdd83e33e8d4811,
  0xfdd82273e8d7860,
  0xfdd806a7e8da8af,
  0xfdd7eacfe8dd8fc,
  0xfdd7cef3e8e0947,
  0xfdd7b30be8e3991,
  0xfdd79717e8e69da,
  0xfdd77b1be8e9a22,
  0xfdd75f17e8eca68,
  0xfdd7430be8efaad,
  0xfdd726f3e8f2af0,
  0xfdd70acfe8f5b32,
  0xfdd6eea7e8f8b73,
  0xfdd6d273e8fbbb2,
  0xfdd6b633e8febf0,
  0xfdd699efe901c2c,
  0xfdd67d9be904c67,
  0xfdd66143e907ca1,
  0xfdd644dfe90acd9,
  0xfdd62873e90dd10,
  0xfdd60bfbe910d46,
  0xfdd5ef7be913d7a,
  0xfdd5d2f3e916dac,
  0xfdd5b65fe919dde,
  0xfdd599c3e91ce0e,
  0xfdd57d1fe91fe3b,
  0xfdd5606fe922e68,
  0xfdd543b7e925e94,
  0xfdd526f3e928ebe,
  0xfdd50a2be92bee7,
  0xfdd4ed53e92ef0f,
  0xfdd4d077e931f35,
  0xfdd4b38fe934f59,
  0xfdd4969be937f7c,
  0xfdd479a3e93af9e,
  0xfdd45c9fe93dfbe,
  0xfdd43f8fe940fdd,
  0xfdd42277e943ffb,
  0xfdd40557e947017,
  0xfdd3e82fe94a031,
  0xfdd3cafbe94d04b,
  0xfdd3adbfe950062,
  0xfdd39077e953079,
  0xfdd37327e95608e,
  0xfdd355cfe9590a1,
  0xfdd3386be95c0b3,
  0xfdd31affe95f0c4,
  0xfdd2fd8be9620d3,
  0xfdd2e00be9650e0,
  0xfdd2c283e9680ed,
  0xfdd2a4efe96b0f7,
  0xfdd28753e96e101,
  0xfdd269afe971109,
  0xfdd24c03e97410e,
  0xfdd22e4be977113,
  0xfdd2108be97a116,
  0xfdd1f2bfe97d118,
  0xfdd1d4ebe980119,
  0xfdd1b70fe983118,
  0xfdd19927e986116,
  0xfdd17b37e989112,
  0xfdd15d3fe98c10d,
  0xfdd13f3be98f106,
  0xfdd1212fe9920fe,
  0xfdd10317e9950f4,
  0xfdd0e4fbe9980e9,
  0xfdd0c6cfe99b0dc,
  0xfdd0a89fe99e0ce,
  0xfdd08a63e9a10be,
  0xfdd06c1fe9a40ad,
  0xfdd04dcfe9a709a,
  0xfdd02f77e9aa086,
  0xfdd01117e9ad071,
  0xfdcff2abe9b005a,
  0xfdcfd437e9b3041,
  0xfdcfb5bbe9b6027,
  0xfdcf9733e9b900b,
  0xfdcf78a3e9bbfee,
  0xfdcf5a0be9befd0,
  0xfdcf3b67e9c1fb0,
  0xfdcf1cbbe9c4f8d,
  0xfdcefe07e9c7f6a,
  0xfdcedf47e9caf45,
  0xfdcec07fe9cdf1f,
  0xfdcea1abe9d0ef8,
  0xfdce82d3e9d3ece,
  0xfdce63ebe9d6ea4,
  0xfdce44ffe9d9e78,
  0xfdce2607e9dce4a,
  0xfdce0707e9dfe1b,
  0xfdcde7fbe9e2dea,
  0xfdcdc8e7e9e5db8,
  0xfdcda9cbe9e8d84,
  0xfdcd8aa3e9ebd4f,
  0xfdcd6b73e9eed18,
  0xfdcd4c3be9f1cdf,
  0xfdcd2cfbe9f4ca5,
  0xfdcd0dafe9f7c6a,
  0xfdccee57e9fac2d,
  0xfdcccef7e9fdbee,
  0xfdccaf8fea00bae,
  0xfdcc901fea03b6d,
  0xfdcc70a3ea06b2a,
  0xfdcc511fea09ae5,
  0xfdcc3193ea0ca9f,
  0xfdcc11fbea0fa57,
  0xfdcbf25bea12a0e,
  0xfdcbd2b3ea159c2,
  0xfdcbb2ffea18975,
  0xfdcb9343ea1b927,
  0xfdcb737fea1e8d8,
  0xfdcb53afea21887,
  0xfdcb33d7ea24834,
  0xfdcb13f3ea277e0,
  0xfdcaf407ea2a78a,
  0xfdcad413ea2d732,
  0xfdcab417ea306d9,
  0xfdca940fea3367f,
  0xfdca73ffea36623,
  0xfdca53e3ea395c5,
  0xfdca33bfea3c566,
  0xfdca1393ea3f505,
  0xfdc9f35fea424a3,
  0xfdc9d31fea4543f,
  0xfdc9b2d7ea483d9,
  0xfdc99283ea4b372,
  0xfdc97227ea4e309,
  0xfdc951c3ea5129f,
  0xfdc93153ea54233,
  0xfdc910dfea571c5,
  0xfdc8f05bea5a156,
  0xfdc8cfd3ea5d0e5,
  0xfdc8af3fea60073,
  0xfdc88ea3ea62fff,
  0xfdc86dfbea65f8a,
  0xfdc84d4bea68f12,
  0xfdc82c93ea6be99,
  0xfdc80bd3ea6ee1f,
  0xfdc7eb07ea71da3,
  0xfdc7ca33ea74d25,
  0xfdc7a953ea77ca6,
  0xfdc7886bea7ac25,
  0xfdc7677bea7dba3,
  0xfdc7467fea80b1f,
  0xfdc7257fea83a99,
  0xfdc7046fea86a12,
  0xfdc6e35bea89989,
  0xfdc6c23bea8c8fe,
  0xfdc6a113ea8f872,
  0xfdc67fdfea927e5,
  0xfdc65ea3ea95755,
  0xfdc63d5fea986c4,
  0xfdc61c13ea9b631,
  0xfdc5fabbea9e59d,
  0xfdc5d95beaa1507,
  0xfdc5b7efeaa4470,
  0xfdc5967beaa73d6,
  0xfdc574ffeaaa33c,
  0xfdc5537beaad29f,
  0xfdc531ebeab0201,
  0xfdc51053eab3161,
  0xfdc4eeb3eab60c0,
  0xfdc4cd07eab901c,
  0xfdc4ab53eabbf77,
  0xfdc48993eabeed1,
  0xfdc467cfeac1e28,
  0xfdc445ffeac4d7f,
  0xfdc42423eac7cd3,
  0xfdc4023feacac26,
  0xfdc3e053eacdb78,
  0xfdc3be5fead0ac7,
  0xfdc39c5fead3a15,
  0xfdc37a57ead6962,
  0xfdc35847ead98ac,
  0xfdc3362beadc7f5,
  0xfdc31407eadf73c,
  0xfdc2f1dbeae2682,
  0xfdc2cfa7eae55c6,
  0xfdc2ad67eae8508,
  0xfdc28b1beaeb449,
  0xfdc268cbeaee388,
  0xfdc2466feaf12c5,
  0xfdc2240beaf4200,
  0xfdc2019beaf713a,
  0xfdc1df23eafa072,
  0xfdc1bca3eafcfa9,
  0xfdc19a1beaffedd,
  0xfdc17787eb02e10,
  0xfdc154ebeb05d42,
  0xfdc13243eb08c71,
  0xfdc10f97eb0bb9e,
  0xfdc0ecdfeb0eaca,
  0xfdc0ca1beb119f5,
  0xfdc0a753eb1491e,
  0xfdc0847feb17845,
  0xfdc0619feb1a76a,
  0xfdc03ebbeb1d68e,
  0xfdc01bcbeb205b0,
  0xfdbff8d3eb234d0,
  0xfdbfd5cfeb263ef,
  0xfdbfb2c3eb2930c,
  0xfdbf8fafeb2c227,
  0xfdbf6c93eb2f140,
  0xfdbf496beb32058,
  0xfdbf263beb34f6e,
  0xfdbf02ffeb37e82,
  0xfdbedfbbeb3ad94,
  0xfdbebc6feb3dca5,
  0xfdbe991beb40bb4,
  0xfdbe75bbeb43ac1,
  0xfdbe5253eb469cd,
  0xfdbe2ee3eb498d7,
  0xfdbe0b67eb4c7df,
  0xfdbde7e3eb4f6e5,
  0xfdbdc457eb525ea,
  0xfdbda0c3eb554ec,
  0xfdbd7d23eb583ed,
  0xfdbd597beb5b2ec,
  0xfdbd35c7eb5e1e9,
  0xfdbd120feb610e5,
  0xfdbcee4beb63fdf,
  0xfdbcca7beb66ed8,
  0xfdbca6a7eb69dce,
  0xfdbc82c7eb6ccc3,
  0xfdbc5edbeb6fbb6,
  0xfdbc3aebeb72aa7,
  0xfdbc16efeb75997,
  0xfdbbf2ebeb78884,
  0xfdbbcedbeb7b770,
  0xfdbbaac3eb7e65b,
  0xfdbb86a3eb81543,
  0xfdbb627beb8442a,
  0xfdbb3e47eb8730e,
  0xfdbb1a0beb8a1f2,
  0xfdbaf5c7eb8d0d3,
  0xfdbad177eb8ffb2,
  0xfdbaad1feb92e90,
  0xfdba88bfeb95d6c,
  0xfdba6453eb98c46,
  0xfdba3fe3eb9bb1f,
  0xfdba1b63eb9e9f5,
  0xfdb9f6dfeba18ca,
  0xfdb9d24feba479d,
  0xfdb9adb7eba766e,
  0xfdb98917ebaa53d,
  0xfdb9646febad40a,
  0xfdb93fbbebb02d6,
  0xfdb91afbebb31a0,
  0xfdb8f637ebb6068,
  0xfdb8d167ebb8f2f,
  0xfdb8ac8febbbdf3,
  0xfdb887afebbecb6,
  0xfdb862c3ebc1b77,
  0xfdb83dcfebc4a36,
  0xfdb818d3ebc78f4,
  0xfdb7f3cbebca7af,
  0xfdb7cebbebcd669,
  0xfdb7a9a3ebd0521,
  0xfdb78483ebd33d7,
  0xfdb75f57ebd628b,
  0xfdb73a23ebd913d,
  0xfdb714e7ebdbfee,
  0xfdb6ef9febdee9d,
  0xfdb6ca4febe1d4a,
  0xfdb6a4f7ebe4bf5,
  0xfdb67f93ebe7a9e,
  0xfdb65a2bebea946,
  0xfdb634b3ebed7eb,
  0xfdb60f37ebf068f,
  0xfdb5e9afebf3531,
  0xfdb5c423ebf63d1,
  0xfdb59e87ebf926f,
  0xfdb578e7ebfc10b,
  0xfdb5533bebfefa5,
  0xfdb52d87ec01e3e,
  0xfdb507cbec04cd5,
  0xfdb4e203ec07b6a,
  0xfdb4bc33ec0a9fd,
  0xfdb4965bec0d88e,
  0xfdb47077ec1071d,
  0xfdb44a8fec135ab,
  0xfdb42497ec16437,
  0xfdb3fe9bec192c0,
  0xfdb3d893ec1c148,
  0xfdb3b283ec1efce,
  0xfdb38c6bec21e53,
  0xfdb3664bec24cd5,
  0xfdb3401fec27b55,
  0xfdb319ebec2a9d4,
  0xfdb2f3abec2d851,
  0xfdb2cd67ec306cc,
  0xfdb2a717ec33545,
  0xfdb280bfec363bc,
  0xfdb25a5bec39231,
  0xfdb233efec3c0a4,
  0xfdb20d7bec3ef16,
  0xfdb1e6ffec41d85,
  0xfdb1c077ec44bf3,
  0xfdb199ebec47a5e,
  0xfdb17353ec4a8c7,
  0xfdb14cafec4d72f,
  0xfdb12603ec50595,
  0xfdb0ff4fec533fa,
  0xfdb0d893ec5625c,
  0xfdb0b1cfec590bc,
  0xfdb08affec5bf1b,
  0xfdb06427ec5ed77,
  0xfdb03d43ec61bd2,
  0xfdb0165bec64a2b,
  0xfdafef67ec67882,
  0xfdafc86bec6a6d6,
  0xfdafa163ec6d529,
  0xfdaf7a53ec7037b,
  0xfdaf533bec731ca,
  0xfdaf2c1bec76017,
  0xfdaf04efec78e62,
  0xfdaeddbfec7bcac,
  0xfdaeb683ec7eaf3,
  0xfdae8f3bec81939,
  0xfdae67ebec8477c,
  0xfdae4097ec875be,
  0xfdae1933ec8a3fe,
  0xfdadf1cbec8d23c,
  0xfdadca57ec90077,
  0xfdada2dbec92eb1,
  0xfdad7b57ec95ce9,
  0xfdad53c7ec98b1f,
  0xfdad2c33ec9b953,
  0xfdad0493ec9e785,
  0xfdacdce7eca15b5,
  0xfdacb537eca43e3,
  0xfdac8d7beca7210,
  0xfdac65b7ecaa03a,
  0xfdac3de7ecace63,
  0xfdac1613ecafc89,
  0xfdabee33ecb2aae,
  0xfdabc647ecb58d0,
  0xfdab9e57ecb86f1,
  0xfdab765becbb50f,
  0xfdab4e57ecbe32c,
  0xfdab264becc1147,
  0xfdaafe33ecc3f60,
  0xfdaad617ecc6d76,
  0xfdaaadefecc9b8b,
  0xfdaa85bbeccc99e,
  0xfdaa5d83eccf7af,
  0xfdaa353fecd25be,
  0xfdaa0cf3ecd53cb,
  0xfda9e49becd81d6,
  0xfda9bc3fecdafdf,
  0xfda993d7ecddde6,
  0xfda96b67ece0beb,
  0xfda942ebece39ee,
  0xfda91a6bece67ef,
  0xfda8f1dfece95ed,
  0xfda8c94becec3ea,
  0xfda8a0abecef1e5,
  0xfda87807ecf1fde,
  0xfda84f57ecf4dd5,
  0xfda8269fecf7bca,
  0xfda7fddbecfa9bd,
  0xfda7d50fecfd7ae,
  0xfda7ac3bed0059d,
  0xfda7835fed0338a,
  0xfda75a7bed06175,
  0xfda7318bed08f5e,
  0xfda70893ed0bd45,
  0xfda6df93ed0eb2a,
  0xfda6b687ed1190d,
  0xfda68d77ed146ee,
  0xfda6645bed174cd,
  0xfda63b33ed1a2aa,
  0xfda61207ed1d085,
  0xfda5e8cfed1fe5e,
  0xfda5bf8fed22c35,
  0xfda59647ed25a0a,
  0xfda56cf3ed287dc,
  0xfda5439bed2b5ad,
  0xfda51a37ed2e37c,
  0xfda4f0cbed31149,
  0xfda4c753ed33f13,
  0xfda49dd7ed36cdb,
  0xfda4744fed39aa2,
  0xfda44abbed3c866,
  0xfda42123ed3f629,
  0xfda3f77fed423e9,
  0xfda3cdd3ed451a8,
  0xfda3a41fed47f64,
  0xfda37a5fed4ad1e,
  0xfda3509bed4dad7,
  0xfda326cbed5088d,
  0xfda2fcf3ed53641,
  0xfda2d30fed563f3,
  0xfda2a927ed591a3,
  0xfda27f33ed5bf52,
  0xfda25533ed5ecfd,
  0xfda22b2fed61aa7,
  0xfda2011fed6484f,
  0xfda1d70bed675f5,
  0xfda1ace7ed6a399,
  0xfda182bfed6d13a,
  0xfda1588bed6feda,
  0xfda12e53ed72c77,
  0xfda1040fed75a13,
  0xfda0d9bfed787ac,
  0xfda0af6bed7b543,
  0xfda0850bed7e2d9,
  0xfda05aa3ed8106c,
  0xfda03033ed83dfd,
  0xfda005b7ed86b8b,
  0xfd9fdb37ed89918,
  0xfd9fb0abed8c6a3,
  0xfd9f8613ed8f42b,
  0xfd9f5b77ed921b2,
  0xfd9f30cfed94f36,
  0xfd9f061fed97cb9,
  0xfd9edb67ed9aa39,
  0xfd9eb0a7ed9d7b7,
  0xfd9e85dbeda0533,
  0xfd9e5b07eda32ad,
  0xfd9e302beda6025,
  0xfd9e0547eda8d9b,
  0xfd9dda57edabb0f,
  0xfd9daf63edae880,
  0xfd9d8463edb15f0,
  0xfd9d5957edb435d,
  0xfd9d2e47edb70c8,
  0xfd9d032bedb9e32,
  0xfd9cd807edbcb99,
  0xfd9cacdbedbf8fe,
  0xfd9c81a7edc2660,
  0xfd9c5667edc53c1,
  0xfd9c2b1fedc811f,
  0xfd9bffcfedcae7c,
  0xfd9bd477edcdbd6,
  0xfd9ba913edd092e,
  0xfd9b7dabedd3683,
  0xfd9b5237edd63d7,
  0xfd9b26bbedd9129,
  0xfd9afb33eddbe79,
  0xfd9acfa7eddebc6,
  0xfd9aa40fede1912,
  0xfd9a786fede465b,
  0xfd9a4cc3ede73a2,
  0xfd9a2113edea0e7,
  0xfd99f557edece2a,
  0xfd99c993edefb6a,
  0xfd999dc7edf28a9,
  0xfd9971f3edf55e5,
  0xfd994613edf831f,
  0xfd991a2bedfb057,
  0xfd98ee3bedfdd8d,
  0xfd98c243ee00ac1,
  0xfd98963fee037f2,
  0xfd986a37ee06522,
  0xfd983e23ee0924f,
  0xfd981203ee0bf7a,
  0xfd97e5dfee0eca3,
  0xfd97b9afee119ca,
  0xfd978d7bee146ee,
  0xfd97613bee17411,
  0xfd9734efee1a131,
  0xfd97089fee1ce4f,
  0xfd96dc43ee1fb6b,
  0xfd96afe3ee22883,
  0xfd968373ee2559b,
  0xfd9656ffee282b0,
  0xfd962a83ee2afc3,
  0xfd95fdfbee2dcd4,
  0xfd95d16bee309e3,
  0xfd95a4d3ee336f0,
  0xfd95782fee363fa,
  0xfd954b87ee39102,
  0xfd951ed3ee3be08,
  0xfd94f217ee3eb0c,
  0xfd94c553ee4180e,
  0xfd949883ee4450d,
  0xfd946bafee4720a,
  0xfd943ecfee49f05,
  0xfd9411e7ee4cbfe,
  0xfd93e4f7ee4f8f5,
  0xfd93b7fbee525e9,
  0xfd938af7ee552db,
  0xfd935defee57fcb,
  0xfd9330d7ee5acb9,
  0xfd9303bbee5d9a5,
  0xfd92d697ee6068e,
  0xfd92a967ee63375,
  0xfd927c2fee6605a,
  0xfd924eefee68d3d,
  0xfd9221a7ee6ba1d,
  0xfd91f453ee6e6fb,
  0xfd91c6fbee713d7,
  0xfd919997ee740b0,
  0xfd916c2bee76d88,
  0xfd913eb3ee79a5d,
  0xfd911137ee7c730,
  0xfd90e3afee7f401,
  0xfd90b61fee820d0,
  0xfd908887ee84d9c,
  0xfd905ae7ee87a66,
  0xfd902d3bee8a72e,
  0xfd8fff87ee8d3f4,
  0xfd8fd1cfee900b7,
  0xfd8fa407ee92d78,
  0xfd8f763bee95a37,
  0xfd8f4867ee986f4,
  0xfd8f1a87ee9b3ae,
  0xfd8eec9fee9e067,
  0xfd8ebeafeea0d1d,
  0xfd8e90b7eea39d0,
  0xfd8e62b3eea6682,
  0xfd8e34abeea9331,
  0xfd8e0697eeabfde,
  0xfd8dd87beeaec88,
  0xfd8daa53eeb1931,
  0xfd8d7c27eeb45d7,
  0xfd8d4defeeb727b,
  0xfd8d1fb3eeb9f1b,
  0xfd8cf16beebcbba,
  0xfd8cc31beebf857,
  0xfd8c94bfeec24f2,
  0xfd8c665feec518b,
  0xfd8c37f3eec7e21,
  0xfd8c097feecaab5,
  0xfd8bdb03eecd746,
  0xfd8bac7beed03d6,
  0xfd8b7defeed3063,
  0xfd8b4f57eed5cee,
  0xfd8b20b7eed8976,
  0xfd8af20feedb5fc,
  0xfd8ac35feede280,
  0xfd8a94a3eee0f02,
  0xfd8a65e3eee3b81,
  0xfd8a3717eee67fe,
  0xfd8a0843eee9479,
  0xfd89d967eeec0f1,
  0xfd89aa7feeeed67,
  0xfd897b93eef19db,
  0xfd894c9beef464d,
  0xfd891d9beef72bc,
  0xfd88ee93eef9f29,
  0xfd88bf83eefcb93,
  0xfd889067eeff7fc,
  0xfd886147ef02462,
  0xfd88321bef050c5,
  0xfd8802e7ef07d26,
  0xfd87d3abef0a985,
  0xfd87a467ef0d5e1,
  0xfd877517ef1023c,
  0xfd8745bfef12e94,
  0xfd871663ef15ae9,
  0xfd86e6fbef1873d,
  0xfd86b787ef1b38e,
  0xfd86880fef1dfdd,
  0xfd86588bef20c29,
  0xfd862903ef23873,
  0xfd85f96fef264bb,
  0xfd85c9d3ef29100,
  0xfd859a2bef2bd43,
  0xfd856a7fef2e984,
  0xfd853ac7ef315c3,
  0xfd850b07ef341ff,
  0xfd84db3fef36e38,
  0xfd84ab6fef39a70,
  0xfd847b97ef3c6a5,
  0xfd844bb7ef3f2d7,
  0xfd841bcbef41f08,
  0xfd83ebd7ef44b36,
  0xfd83bbdbef47761,
  0xfd838bd7ef4a38a,
  0xfd835bcbef4cfb1,
  0xfd832bb3ef4fbd6,
  0xfd82fb97ef527f7,
  0xfd82cb6fef55417,
  0xfd829b3fef58034,
  0xfd826b07ef5ac4f,
  0xfd823ac7ef5d868,
  0xfd820a7bef6047e,
  0xfd81da27ef63092,
  0xfd81a9cfef65ca4,
  0xfd81796bef688b3,
  0xfd8148ffef6b4c0,
  0xfd811887ef6e0ca,
  0xfd80e80bef70cd3,
  0xfd80b783ef738d8,
  0xfd8086f3ef764dc,
  0xfd80565bef790dd,
  0xfd8025bbef7bcdb,
  0xfd7ff513ef7e8d7,
  0xfd7fc463ef814d1,
  0xfd7f93a7ef840c8,
  0xfd7f62e3ef86cbd,
  0xfd7f321bef898b0,
  0xfd7f0143ef8c4a0,
  0xfd7ed067ef8f08e,
  0xfd7e9f83ef91c79,
  0xfd7e6e93ef94862,
  0xfd7e3d9fef97449,
  0xfd7e0c9fef9a02d,
  0xfd7ddb97ef9cc0f,
  0xfd7daa87ef9f7ee,
  0xfd7d796fefa23cb,
  0xfd7d484befa4fa6,
  0xfd7d1723efa7b7e,
  0xfd7ce5efefaa754,
  0xfd7cb4b3efad327,
  0xfd7c836fefafef8,
  0xfd7c5223efb2ac6,
  0xfd7c20cbefb5692,
  0xfd7bef6fefb825c,
  0xfd7bbe07efbae23,
  0xfd7b8c97efbd9e8,
  0xfd7b5b1fefc05aa,
  0xfd7b29a3efc3169,
  0xfd7af81befc5d26,
  0xfd7ac687efc88e1,
  0xfd7a94efefcb49a,
  0xfd7a634befce050,
  0xfd7a319fefd0c04,
  0xfd79ffebefd37b5,
  0xfd79ce2fefd6364,
  0xfd799c6befd8f10,
  0xfd796a9fefdbaba,
  0xfd7938c7efde662,
  0xfd7906e7efe1207,
  0xfd78d4ffefe3daa,
  0xfd78a313efe694a,
  0xfd787117efe94e8,
  0xfd783f17efec083,
  0xfd780d0fefeec1c,
  0xfd77dafbeff17b2,
  0xfd77a8e3eff4346,
  0xfd7776bfeff6ed7,
  0xfd774493eff9a66,
  0xfd77125feffc5f3,
  0xfd76e023efff17d,
  0xfd76addbf000e82,
  0xfd767b8ff002445,
  0xfd764937f003a06,
  0xfd7616d7f004fc6,
  0xfd75e46ff006585,
  0xfd75b1fff007b43,
  0xfd757f87f0090ff,
  0xfd754d07f00a6ba,
  0xfd751a7bf00bc74,
  0xfd74e7ebf00d22d,
  0xfd74b54ff00e7e4,
  0xfd7482abf00fd9a,
  0xfd744ffff01134f,
  0xfd741d4bf012903,
  0xfd73ea8ff013eb6,
  0xfd73b7cbf015467,
  0xfd7384fbf016a17,
  0xfd735227f017fc5,
  0xfd731f47f019573,
  0xfd72ec5ff01ab1f,
  0xfd72b96ff01c0ca,
  0xfd728677f01d674,
  0xfd725377f01ec1c,
  0xfd72206bf0201c3,
  0xfd71ed5bf021769,
  0xfd71ba3ff022d0e,
  0xfd71871bf0242b2,
  0xfd7153eff025854,
  0xfd7120bff026df5,
  0xfd70ed7ff028394,
  0xfd70ba3bf029933,
  0xfd7086eff02aed0,
  0xfd70539bf02c46b,
  0xfd70203bf02da06,
  0xfd6fecd3f02ef9f,
  0xfd6fb967f030537,
  0xfd6f85eff031ace,
  0xfd6f526ff033064,
  0xfd6f1ee3f0345f8,
  0xfd6eeb53f035b8b,
  0xfd6eb7bbf03711d,
  0xfd6e8417f0386ad,
  0xfd6e506bf039c3d,
  0xfd6e1cbbf03b1cb,
  0xfd6de8fff03c757,
  0xfd6db53bf03dce3,
  0xfd6d816bf03f26d,
  0xfd6d4d97f0407f6,
  0xfd6d19bbf041d7d,
  0xfd6ce5d3f043304,
  0xfd6cb1e7f044889,
  0xfd6c7deff045e0d,
  0xfd6c49eff04738f,
  0xfd6c15e7f048910,
  0xfd6be1d7f049e90,
  0xfd6badbff04b40f,
  0xfd6b799ff04c98d,
  0xfd6b4573f04df09,
  0xfd6b1143f04f484,
  0xfd6add07f0509fd,
  0xfd6aa8c7f051f76,
  0xfd6a747bf0534ed,
  0xfd6a4027f054a62,
  0xfd6a0bcbf055fd7,
  0xfd69d767f05754a,
  0xfd69a2fbf058abc,
  0xfd696e83f05a02c,
  0xfd693a07f05b59c,
  0xfd69057ff05cb0a,
  0xfd68d0f3f05e076,
  0xfd689c5bf05f5e2,
  0xfd6867bbf060b4c,
  0xfd683313f0620b5,
  0xfd67fe63f06361c,
  0xfd67c9abf064b83,
  0xfd6794ebf0660e8,
  0xfd67601ff06764b,
  0xfd672b4ff068bae,
  0xfd66f673f06a10f,
  0xfd66c193f06b66f,
  0xfd668ca7f06cbcd,
  0xfd6657b3f06e12a,
  0xfd6622b7f06f686,
  0xfd65edb3f070be1,
  0xfd65b8a7f07213a,
  0xfd658393f073692,
  0xfd654e73f074be9,
  0xfd65194ff07613d,
  0xfd64e423f077691,
  0xfd64aeebf078be4,
  0xfd6479abf07a135,
  0xfd644467f07b686,
  0xfd640f17f07cbd4,
  0xfd63d9bff07e122,
  0xfd63a45ff07f66e,
  0xfd636ef7f080bb9,
  0xfd633983f082103,
  0xfd63040bf08364b,
  0xfd62ce87f084b92,
  0xfd6298fff0860d8,
  0xfd62636bf08761c,
  0xfd622dd3f088b5f,
  0xfd61f82ff08a0a1,
  0xfd61c283f08b5e1,
  0xfd618ccff08cb20,
  0xfd615713f08e05e,
  0xfd61214ff08f59b,
  0xfd60eb83f090ad6,
  0xfd60b5abf09200f,
  0xfd607fcff093548,
  0xfd6049e7f094a7f,
  0xfd6013fbf095fb5,
  0xfd5fde03f0974e9,
  0xfd5fa807f098a1d,
  0xfd5f71fff099f4f,
  0xfd5f3beff09b47f,
  0xfd5f05d7f09c9ae,
  0xfd5ecfb7f09dedc,
  0xfd5e998ff09f409,
  0xfd5e635ff0a0934,
  0xfd5e2d23f0a1e5e,
  0xfd5df6e3f0a3386,
  0xfd5dc09bf0a48ae,
  0xfd5d8a47f0a5dd3,
  0xfd5d53ebf0a72f8,
  0xfd5d1d8bf0a881b,
  0xfd5ce71ff0a9d3d,
  0xfd5cb0abf0ab25e,
  0xfd5c7a2ff0ac77d,
  0xfd5c43abf0adc9b,
  0xfd5c0d1ff0af1b7,
  0xfd5bd68bf0b06d3,
  0xfd5b9feff0b1bec,
  0xfd5b694bf0b3105,
  0xfd5b329ff0b461c,
  0xfd5afbe7f0b5b32,
  0xfd5ac52bf0b7047,
  0xfd5a8e63f0b855a,
  0xfd5a5797f0b9a6b,
  0xfd5a20bff0baf7c,
  0xfd59e9dff0bc48b,
  0xfd59b2fbf0bd999,
  0xfd597c0bf0beea4,
  0xfd594513f0c03b0,
  0xfd590e13f0c18b9,
  0xfd58d70bf0c2dc2,
  0xfd589ffbf0c42c9,
  0xfd5868e3f0c57ce,
  0xfd5831c3f0c6cd3,
  0xfd57fa97f0c81d6,
  0xfd57c367f0c96d7,
  0xfd578c2ff0cabd7,
  0xfd5754ebf0cc0d6,
  0xfd571da3f0cd5d4,
  0xfd56e64ff0cead0,
  0xfd56aef3f0cffcb,
  0xfd567793f0d14c4,
  0xfd564027f0d29bc,
  0xfd5608b3f0d3eb3,
  0xfd55d137f0d53a9,
  0xfd5599b3f0d689d,
  0xfd556227f0d7d8f,
  0xfd552a93f0d9281,
  0xfd54f2f7f0da771,
  0xfd54bb53f0dbc5f,
  0xfd5483a3f0dd14c,
  0xfd544beff0de638,
  0xfd541433f0dfb23,
  0xfd53dc6bf0e100c,
  0xfd53a49ff0e24f3,
  0xfd536cc7f0e39da,
  0xfd5334ebf0e4ebf,
  0xfd52fd03f0e63a2,
  0xfd52c513f0e7885,
  0xfd528d1ff0e8d66,
  0xfd52551ff0ea245,
  0xfd521d17f0eb723,
  0xfd51e507f0ecc00,
  0xfd51aceff0ee0db,
  0xfd5174cff0ef5b5,
  0xfd513ca7f0f0a8e,
  0xfd510477f0f1f65,
  0xfd50cc3ff0f343b,
  0xfd5093fff0f490f,
  0xfd505bb7f0f5de3,
  0xfd502367f0f72b4,
  0xfd4feb0bf0f8785,
  0xfd4fb2abf0f9c54,
  0xfd4f7a43f0fb121,
  0xfd4f41cff0fc5ed,
  0xfd4f0957f0fdab8,
  0xfd4ed0d3f0fef82,
  0xfd4e984bf10044a,
  0xfd4e5fb7f101910,
  0xfd4e271ff102dd5,
  0xfd4dee7bf104299,
  0xfd4db5cff10575c,
  0xfd4d7d1ff106c1c,
  0xfd4d4463f1080dc,
  0xfd4d0ba3f10959a,
  0xfd4cd2d7f10aa57,
  0xfd4c9a03f10bf12,
  0xfd4c6127f10d3cc,
  0xfd4c2843f10e885,
  0xfd4bef57f10fd3c,
  0xfd4bb663f1111f2,
  0xfd4b7d67f1126a7,
  0xfd4b4463f113b5a,
  0xfd4b0b57f11500c,
  0xfd4ad243f1164bc,
  0xfd4a9927f11796b,
  0xfd4a5ffff118e19,
  0xfd4a26d3f11a2c5,
  0xfd49ed9ff11b76f,
  0xfd49b463f11cc19,
  0xfd497b1bf11e0c1,
  0xfd4941cff11f567,
  0xfd49087bf120a0c,
  0xfd48cf1bf121eb0,
  0xfd4895b7f123352,
  0xfd485c47f1247f3,
  0xfd4822d3f125c92,
  0xfd47e953f127130,
  0xfd47afcff1285cd,
  0xfd47763ff129a68,
  0xfd473cabf12af02,
  0xfd47030bf12c39a,
  0xfd46c967f12d831,
  0xfd468fb7f12ecc7,
  0xfd4655fff13015b,
  0xfd461c43f1315ee,
  0xfd45e27bf132a7f,
  0xfd45a8abf133f0f,
  0xfd456ed7f13539e,
  0xfd4534f7f13682b,
  0xfd44fb0ff137cb6,
  0xfd44c11ff139140,
  0xfd44872bf13a5c9,
  0xfd444d2bf13ba51,
  0xfd441323f13ced6,
  0xfd43d913f13e35b,
  0xfd439efbf13f7de,
  0xfd4364dbf140c60,
  0xfd432ab7f1420e0,
  0xfd42f087f14355f,
  0xfd42b64ff1449dc,
  0xfd427c0ff145e58,
  0xfd4241c7f1472d3,
  0xfd420777f14874c,
  0xfd41cd1ff149bc3,
  0xfd4192bff14b03a,
  0xfd41585bf14c4ae,
  0xfd411debf14d921,
  0xfd40e373f14ed93,
  0xfd40a8f3f150203,
  0xfd406e6bf151672,
  0xfd4033dbf152ae0,
  0xfd3ff943f153f4c,
  0xfd3fbea3f1553b7,
  0xfd3f83fbf156820,
  0xfd3f494bf157c88,
  0xfd3f0e93f1590ef,
  0xfd3ed3d3f15a554,
  0xfd3e990bf15b9b7,
  0xfd3e5e3bf15ce19,
  0xfd3e2363f15e27a,
  0xfd3de883f15f6d9,
  0xfd3dad9bf160b37,
  0xfd3d72abf161f93,
  0xfd3d37b3f1633ee,
  0xfd3cfcb3f164847,
  0xfd3cc1abf165c9f,
  0xfd3c8697f1670f6,
  0xfd3c4b7ff16854b,
  0xfd3c105ff16999f,
  0xfd3bd537f16adf1,
  0xfd3b9a07f16c241,
  0xfd3b5ecff16d691,
  0xfd3b2393f16eade,
  0xfd3ae84bf16ff2b,
  0xfd3aacfbf171376,
  0xfd3a71a3f1727bf,
  0xfd3a3643f173c07,
  0xfd39fadbf17504e,
  0xfd39bf6bf176493,
  0xfd3983f3f1778d6,
  0xfd394873f178d19,
  0xfd390cebf17a159,
  0xfd38d15bf17b598,
  0xfd3895c3f17c9d6,
  0xfd385a23f17de13,
  0xfd381e7ff17f24d,
  0xfd37e2cff180687,
  0xfd37a717f181abf,
  0xfd376b57f182ef5,
  0xfd372f8ff18432a,
  0xfd36f3c3f18575e,
  0xfd36b7ebf186b90,
  0xfd367c0bf187fc0,
  0xfd364023f1893ef,
  0xfd360437f18a81d,
  0xfd35c83ff18bc49,
  0xfd358c3ff18d074,
  0xfd35503bf18e49d,
  0xfd35142bf18f8c5,
  0xfd34d813f190ceb,
  0xfd349bf7f19210f,
  0xfd345fd3f193532,
  0xfd3423a3f194954,
  0xfd33e76ff195d75,
  0xfd33ab2ff197193,
  0xfd336eebf1985b1,
  0xfd33329bf1999cd,
  0xfd32f647f19ade7,
  0xfd32b9ebf19c200,
  0xfd327d83f19d618,
  0xfd324117f19ea2e,
  0xfd3204a3f19fe42,
  0xfd31c823f1a1255,
  0xfd318b9ff1a2667,
  0xfd314f13f1a3a77,
  0xfd31127ff1a4e85,
  0xfd30d5e3f1a6292,
  0xfd30993ff1a769e,
  0xfd305c8ff1a8aa8,
  0xfd301fdbf1a9eb1,
  0xfd2fe31ff1ab2b8,
  0xfd2fa65bf1ac6be,
  0xfd2f698ff1adac2,
  0xfd2f2cbff1aeec4,
  0xfd2eefe3f1b02c5,
  0xfd2eb2fff1b16c5,
  0xfd2e7613f1b2ac3,
  0xfd2e391ff1b3ec0,
  0xfd2dfc27f1b52bb,
  0xfd2dbf23f1b66b5,
  0xfd2d8217f1b7aad,
  0xfd2d4507f1b8ea4,
  0xfd2d07ebf1ba299,
  0xfd2ccacbf1bb68d,
  0xfd2c8d9ff1bca7f,
  0xfd2c506ff1bde6f,
  0xfd2c1333f1bf25f,
  0xfd2bd5f3f1c064c,
  0xfd2b98a7f1c1a39,
  0xfd2b5b57f1c2e23,
  0xfd2b1dfff1c420c,
  0xfd2ae09ff1c55f4,
  0xfd2aa337f1c69da,
  0xfd2a65c7f1c7dbf,
  0xfd2a284ff1c91a2,
  0xfd29eacff1ca584,
  0xfd29ad47f1cb964,
  0xfd296fb7f1ccd42,
  0xfd29321ff1ce11f,
  0xfd28f47ff1cf4fb,
  0xfd28b6d7f1d08d5,
  0xfd28792bf1d1cae,
  0xfd283b73f1d3085,
  0xfd27fdb3f1d445a,
  0xfd27bfeff1d582e,
  0xfd278223f1d6c00,
  0xfd27444ff1d7fd1,
  0xfd27066ff1d93a0,
  0xfd26c88bf1da76e,
  0xfd268a9ff1dbb3b,
  0xfd264ca7f1dcf06,
  0xfd260eabf1de2cf,
  0xfd25d0a7f1df697,
  0xfd25929bf1e0a5d,
  0xfd255487f1e1e22,
  0xfd25166bf1e31e5,
  0xfd24d847f1e45a7,
  0xfd249a1ff1e5967,
  0xfd245bebf1e6d26,
  0xfd241daff1e80e3,
  0xfd23df6ff1e949f,
  0xfd23a123f1ea859,
  0xfd2362d3f1ebc12,
  0xfd232477f1ecfc9,
  0xfd22e617f1ee37e,
  0xfd22a7aff1ef732,
  0xfd22693bf1f0ae5,
  0xfd222ac3f1f1e96,
  0xfd21ec43f1f3245,
  0xfd21adbbf1f45f3,
  0xfd216f2bf1f599f,
  0xfd213093f1f6d4a,
  0xfd20f1f3f1f80f3,
  0xfd20b34ff1f949b,
  0xfd20749ff1fa841,
  0xfd2035e7f1fbbe6,
  0xfd1ff72bf1fcf89,
  0xfd1fb863f1fe32b,
  0xfd1f7997f1ff6cb,
  0xfd1f3ac3f200a69,
  0xfd1efbe7f201e06,
  0xfd1ebcfff2031a2,
  0xfd1e7e13f20453b,
  0xfd1e3f1ff2058d4,
  0xfd1e0023f206c6a,
  0xfd1dc11ff208000,
  0xfd1d8217f209393,
  0xfd1d4303f20a726,
  0xfd1d03e7f20bab6,
  0xfd1cc4c7f20ce45,
  0xfd1c859bf20e1d3,
  0xfd1c466bf20f55f,
  0xfd1c0733f2108e9,
  0xfd1bc7eff211c72,
  0xfd1b88a7f212ff9,
  0xfd1b4957f21437f,
  0xfd1b09fff215703,
  0xfd1aca9ff216a86,
  0xfd1a8b3bf217e07,
  0xfd1a4bcbf219186,
  0xfd1a0c57f21a503,
  0xfd19ccd7f21b880,
  0xfd198d53f21cbfb,
  0xfd194dc7f21df74,
  0xfd190e2ff21f2eb,
  0xfd18ce93f220662,
  0xfd188eeff2219d6,
  0xfd184f43f222d49,
  0xfd180f8ff2240bb,
  0xfd17cfd3f22542b,
  0xfd179013f226799,
  0xfd175047f227b06,
  0xfd171073f228e71,
  0xfd16d09bf22a1db,
  0xfd1690bbf22b543,
  0xfd1650cff22c8a9,
  0xfd1610dff22dc0e,
  0xfd15d0e7f22ef71,
  0xfd1590e7f2302d3,
  0xfd1550dff231633,
  0xfd1510cff232992,
  0xfd14d0bbf233cef,
  0xfd14909bf23504b,
  0xfd145077f2363a5,
  0xfd141047f2376fd,
  0xfd13d013f238a54,
  0xfd138fd7f239da9,
  0xfd134f93f23b0fc,
  0xfd130f47f23c44f,
  0xfd12cef3f23d79f,
  0xfd128e97f23eaee,
  0xfd124e33f23fe3b,
  0xfd120dcbf241187,
  0xfd11cd57f2424d1,
  0xfd118cdff24381a,
  0xfd114c5ff244b61,
  0xfd110bd7f245ea6,
  0xfd10cb47f2471ea,
  0xfd108aaff24852c,
  0xfd104a0ff24986d,
  0xfd100967f24abac,
  0xfd0fc8bbf24bee9,
  0xfd0f8803f24d225,
  0xfd0f4747f24e55f,
  0xfd0f0683f24f898,
  0xfd0ec5b7f250bcf,
  0xfd0e84e3f251f05,
  0xfd0e4407f253239,
  0xfd0e0323f25456b,
  0xfd0dc237f25589c,
  0xfd0d8147f256bcb,
  0xfd0d404bf257ef8,
  0xfd0cff4bf259224,
  0xfd0cbe43f25a54f,
  0xfd0c7d33f25b878,
  0xfd0c3c1bf25cb9e,
  0xfd0bfafbf25dec4,
  0xfd0bb9d7f25f1e7,
  0xfd0b78a7f26050a,
  0xfd0b3773f26182b,
  0xfd0af633f262b4a,
  0xfd0ab4eff263e67,
  0xfd0a73a3f265183,
  0xfd0a324ff26649e,
  0xfd09f0f3f2677b6,
  0xfd09af8ff268ace,
  0xfd096e27f269de3,
  0xfd092cb3f26b0f7,
  0xfd08eb3bf26c409,
  0xfd08a9bbf26d71a,
  0xfd086833f26ea29,
  0xfd0826a3f26fd37,
  0xfd07e50bf271043,
  0xfd07a36bf27234d,
  0xfd0761c3f273656,
  0xfd072017f27495d,
  0xfd06de63f275c62,
  0xfd069ca7f276f66,
  0xfd065adff278268,
  0xfd061917f279569,
  0xfd05d743f27a868,
  0xfd059567f27bb65,
  0xfd055387f27ce61,
  0xfd05119bf27e15b,
  0xfd04cfabf27f454,
  0xfd048db3f28074b,
  0xfd044bb3f281a40,
  0xfd0409abf282d34,
  0xfd03c79bf284026,
  0xfd038587f285316,
  0xfd03436bf286605,
  0xfd030143f2878f2,
  0xfd02bf17f288bde,
  0xfd027ce3f289ec8,
  0xfd023aa7f28b1b0,
  0xfd01f867f28c497,
  0xfd01b61bf28d77c,
  0xfd0173cbf28ea5f,
  0xfd01316ff28fd41,
  0xfd00ef0ff291021,
  0xfd00aca7f292300,
  0xfd006a3bf2935dd,
  0xfd0027c3f2948b8,
  0xfcffe543f295b92,
  0xfcffa2bff296e6a,
  0xfcff6033f298140,
  0xfcff1d9ff299415,
  0xfcfedb03f29a6e8,
  0xfcfe985ff29b9ba,
  0xfcfe55b7f29cc89,
  0xfcfe1307f29df57,
  0xfcfdd04bf29f223,
  0xfcfd8d8bf2a04ee,
  0xfcfd4ac3f2a17b8,
  0xfcfd07f7f2a2a7f,
  0xfcfcc51ff2a3d45,
  0xfcfc823ff2a500a,
  0xfcfc3f5bf2a62cd,
  0xfcfbfc6ff2a758e,
  0xfcfbb97bf2a884d,
  0xfcfb767ff2a9b0b,
  0xfcfb337bf2aadc7,
  0xfcfaf073f2ac082,
  0xfcfaad5ff2ad33a,
  0xfcfa6a47f2ae5f2,
  0xfcfa2727f2af8a7,
  0xfcf9e3fff2b0b5b,
  0xfcf9a0cff2b1e0d,
  0xfcf95d9bf2b30be,
  0xfcf91a5bf2b436d,
  0xfcf8d717f2b561a,
  0xfcf893cbf2b68c6,
  0xfcf85077f2b7b70,
  0xfcf80d1bf2b8e19,
  0xfcf7c9bbf2ba0bf,
  0xfcf7864ff2bb364,
  0xfcf742dff2bc608,
  0xfcf6ff67f2bd8aa,
  0xfcf6bbe7f2beb4a,
  0xfcf6785ff2bfde8,
  0xfcf634d3f2c1085,
  0xfcf5f13bf2c2320,
  0xfcf5ad9ff2c35ba,
  0xfcf569fbf2c4851,
  0xfcf5264ff2c5ae7,
  0xfcf4e29ff2c6d7c,
  0xfcf49ee3f2c800f,
  0xfcf45b23f2c92a0,
  0xfcf4175bf2ca52f,
  0xfcf3d38bf2cb7bd,
  0xfcf38fb3f2cca4a,
  0xfcf34bd3f2cdcd4,
  0xfcf307eff2cef5d,
  0xfcf2c403f2d01e4,
  0xfcf2800ff2d146a,
  0xfcf23c13f2d26ed,
  0xfcf1f80ff2d3970,
  0xfcf1b403f2d4bf0,
  0xfcf16ff3f2d5e6f,
  0xfcf12bdbf2d70ec,
  0xfcf0e7bbf2d8368,
  0xfcf0a393f2d95e1,
  0xfcf05f67f2da85a,
  0xfcf01b2ff2dbad0,
  0xfcefd6f7f2dcd44,
  0xfcef92b3f2ddfb7,
  0xfcef4e67f2df229,
  0xfcef0a13f2e0498,
  0xfceec5bbf2e1707,
  0xfcee8157f2e2973,
  0xfcee3ceff2e3bde,
  0xfcedf87ff2e4e47,
  0xfcedb40bf2e60ae,
  0xfced6f8bf2e7314,
  0xfced2b07f2e8578,
  0xfcece67bf2e97da,
  0xfceca1e7f2eaa3b,
  0xfcec5d4bf2ebc9a,
  0xfcec18abf2ecef7,
  0xfcebd3fff2ee153,
  0xfceb8f4ff2ef3ad,
  0xfceb4a97f2f0605,
  0xfceb05d7f2f185b,
  0xfceac113f2f2ab0,
  0xfcea7c47f2f3d03,
  0xfcea376ff2f4f55,
  0xfce9f293f2f61a5,
  0xfce9adb3f2f73f3,
  0xfce968c7f2f863f,
  0xfce923d7f2f988a,
  0xfce8dedff2faad3,
  0xfce899dff2fbd1a,
  0xfce854d7f2fcf60,
  0xfce80fcbf2fe1a4,
  0xfce7cab3f2ff3e6,
  0xfce78597f300627,
  0xfce74073f301865,
  0xfce6fb4bf302aa3,
  0xfce6b617f303cde,
  0xfce670dff304f18,
  0xfce62b9ff306150,
  0xfce5e657f307386,
  0xfce5a10bf3085bb,
  0xfce55bb3f3097ee,
  0xfce51657f30aa1f,
  0xfce4d0f3f30bc4f,
  0xfce48b87f30ce7c,
  0xfce44617f30e0a9,
  0xfce4009bf30f2d3,
  0xfce3bb1bf3104fc,
  0xfce37593f311723,
  0xfce33007f312948,
  0xfce2ea6ff313b6c,
  0xfce2a4d3f314d8e,
  0xfce25f2ff315fae,
  0xfce21983f3171cc,
  0xfce1d3d3f3183e9,
  0xfce18e17f319604,
  0xfce14857f31a81d,
  0xfce10293f31ba34,
  0xfce0bcc3f31cc4a,
  0xfce076eff31de5f,
  0xfce03113f31f071,
  0xfcdfeb2ff320282,
  0xfcdfa543f321491,
  0xfcdf5f4ff32269e,
  0xfcdf1957f3238aa,
  0xfcded357f324ab4,
  0xfcde8d4ff325cbc,
  0xfcde473ff326ec2,
  0xfcde012bf3280c7,
  0xfcddbb0ff3292ca,
  0xfcdd74ebf32a4cc,
  0xfcdd2ebff32b6cb,
  0xfcdce88ff32c8c9,
  0xfcdca257f32dac5,
  0xfcdc5c17f32ecc0,
  0xfcdc15cff32feb8,
  0xfcdbcf7ff3310af,
  0xfcdb892bf3322a5,
  0xfcdb42cff333498,
  0xfcdafc6bf33468a,
  0xfcdab603f33587a,
  0xfcda6f8ff336a68,
  0xfcda2917f337c55,
  0xfcd9e297f338e40,
  0xfcd99c13f33a029,
  0xfcd95583f33b210,
  0xfcd90eeff33c3f6,
  0xfcd8c853f33d5da,
  0xfcd881b3f33e7bc,
  0xfcd83b07f33f99d,
  0xfcd7f457f340b7b,
  0xfcd7ad9ff341d58,
  0xfcd766dff342f34,
  0xfcd7201bf34410d,
  0xfcd6d94ff3452e5,
  0xfcd6927bf3464bb,
  0xfcd64b9ff34768f,
  0xfcd604bff348862,
  0xfcd5bdd7f349a33,
  0xfcd576e7f34ac02,
  0xfcd52feff34bdcf,
  0xfcd4e8f3f34cf9b,
  0xfcd4a1eff34e165,
  0xfcd45ae3f34f32d,
};


#endif


float cosisin_phase_correction_factor()
{
#if defined USE_SIN_COS_LOOKUP_TABLE
  return 4 * NR_ENTRIES / M_PI_F;
#else
  return M_1_PI_F;
#endif
}


float2 cosisin_phase_corrected(float arg)
{
#if defined USE_SIN_COS_LOOKUP_TABLE
  typedef unsigned __attribute__((__ap_int(NR_BITS + 3))) arg_int_t;
  typedef unsigned __attribute__((__ap_int(NR_BITS)))     index_t;

  arg_int_t arg_int = convert_int_rte(arg);
  index_t   index   = arg_int /* & (NR_ENTRIES - 1) */;
  uint3_t   octant  = arg_int >> NR_BITS;

  if (BIT(octant, 0))
    index = ~index & (NR_ENTRIES - 1);

  uint60_t uval = cosisin_table[index];
  float2   val  = (float2) (as_float((uint) (uval >> 30)), as_float((uint) (uint30_t) uval));

  if (BIT(octant, 0) ^ BIT(octant, 1))
    val = val.yx;

  if (BIT(octant, 1) ^ BIT(octant, 2))
    val.x = -val.x;

  if (BIT(octant, 2))
    val.y = -val.y;

  return val;
#else
  return (float2) (cospi(arg), sinpi(arg));
#endif
}


float2 cosisin(float arg)
{
  //return (float2) (cos(arg), sin(arg));
  return cosisin_phase_corrected(arg * cosisin_phase_correction_factor());
}



inline float3 load_float3(__global const volatile float ptr[3])
{
  float3 f;

  for (uint2_t i = 0; i < 3; i ++)
    f[i] = ptr[i];

  return f;
}


inline void do_subgrid(
  float8 subgrid[restrict NR_PIXELS],
  float8 pixel_values[restrict NR_PIXELS],
  const float8 visibilities[restrict NR_CHANNELS],
  const float wave_numbers[restrict NR_CHANNELS],
  const float3 lmn[restrict NR_PIXELS],
  const float3 uvw,
  const float phase_offsets[restrict NR_PIXELS],
  unsigned short time)
{
#if 1
  float phase_indices[NR_PIXELS];

#pragma unroll MAX(UNROLL_PIXELS_FACTOR * UNROLL_VISIBILITIES_FACTOR / NR_CHANNELS, 1)
  for (unsigned short pixel = 0; pixel < NR_PIXELS; pixel ++)
    phase_indices[pixel] = __fpga_reg(-dot(__fpga_reg(lmn[pixel]), __fpga_reg(uvw)));

#pragma loop_coalesce
#pragma ii 1
#pragma ivdep array(pixel_values)
#pragma max_interleaving 1
  for (unsigned short chan_major = 0; chan_major < NR_CHANNELS; chan_major += UNROLL_VISIBILITIES_FACTOR) {
#pragma max_interleaving 1
#pragma unroll UNROLL_PIXELS_FACTOR
    for (unsigned short pixel = 0; pixel < NR_PIXELS; pixel ++) {
      float8 pixel_value = __fpga_reg(pixel_values[pixel]);

      if (time == 0 && chan_major == 0)
	pixel_value = 0;

#pragma unroll
      for (unsigned short chan_minor = 0; chan_minor < UNROLL_VISIBILITIES_FACTOR; chan_minor ++) {
	unsigned short chan       = chan_major + chan_minor;
	float          phase      = __fpga_reg(phase_indices[pixel]) * __fpga_reg(wave_numbers[chan]) + __fpga_reg(phase_offsets[pixel]);
	float2         phasor     = cosisin_phase_corrected(phase);
	float8         visibility = __fpga_reg(visibilities[chan]);

	pixel_value.even +=  phasor.x * visibility.even;
	pixel_value.odd  +=  phasor.x * visibility.odd;
	pixel_value.even += -phasor.y * visibility.odd;
	pixel_value.odd  +=  phasor.y * visibility.even;

	pixel_value = __fpga_reg(pixel_value);
      }

      pixel_values[pixel] = pixel_value;

      if (time == NR_TIMESTEPS - 1 && chan_major == NR_CHANNELS - UNROLL_VISIBILITIES_FACTOR)
	subgrid[pixel] = __fpga_reg(pixel_value);
    }
  }
#else
  for (unsigned short pixel = 0; pixel < NR_PIXELS; pixel ++)
    subgrid[pixel] = pixel;
#endif
}


void apply_aterm_and_spheroidal(
  float8 subgrid_out[restrict NR_PIXELS],
  const float8 subgrids_in[restrict NR_GRIDDERS][NR_PIXELS],
  __global const volatile float8 aterms[restrict NR_STATIONS][NR_PIXELS],
  uint2 stationPair,
  const float spheroidal[restrict NR_PIXELS],
  unsigned short gridder)
{
#if 0
  union matrix _aterms[2][NR_PIXELS];

#pragma ii 1
#pragma loop_coalesce
  for (uint2_t s = 0; s < 2; s ++)
    for (unsigned short pixel = 0; pixel < NR_PIXELS; pixel ++)
      _aterms[s][pixel].f8 = aterms[stationPair[s]][pixel];

  for (unsigned short pixel = 0; pixel < NR_PIXELS; pixel ++)
    subgrid_out[pixel] = matmul(matmul(hermitian(_aterms[0][pixel]), subgrids_in[pixel]), _aterms[1][pixel]) * spheroidal[pixel];
#else
  union matrix tmp[NR_PIXELS] __attribute__((private_copies(2)));

#pragma ii 1
#pragma loop_coalesce
  for (unsigned short pixel = 0; pixel < NR_PIXELS; pixel ++) {
    union matrix aterm, pixel_values[NR_GRIDDERS] __attribute__((register));
    aterm.f8 = aterms[stationPair[0]][pixel];

#pragma unroll
    for (unsigned const_gridder = 0; const_gridder < NR_GRIDDERS; const_gridder ++)
      pixel_values[const_gridder].f8 = subgrids_in[const_gridder][pixel];

    for (uint2_t y = 0; y < 2; y ++)
      for (uint2_t x = 0; x < 2; x ++) {
	float2 sum = 0;

#pragma unroll
	for (uint2_t k = 0; k < 2; k ++)
	  sum += cmul(hermitian(aterm).f2[y][k], pixel_values[gridder].f2[k][x]);

	tmp[pixel].f2[y][x] = sum;
      }
  }

#pragma ii 1
#pragma loop_coalesce
  for (unsigned short pixel = 0; pixel < NR_PIXELS; pixel ++) {
    union matrix aterm, out;
    aterm.f8 = aterms[stationPair[1]][pixel];

    for (uint2_t y = 0; y < 2; y ++)
      for (uint2_t x = 0; x < 2; x ++) {
	float2 sum = 0;

#pragma unroll
	for (uint2_t k = 0; k < 2; k ++)
	  sum += cmul(tmp[pixel].f2[y][k], aterm.f2[k][x]);

	out.f2[y][x] = sum * spheroidal[pixel];

	if (y == 1 && x == 1)
	  subgrid_out[pixel] = out.f8;
      }
  }
#endif
}


#define FFT_FORWARD -1
#define FFT_BACKWARD 1



__constant float2 cosisin_16[16] =
{ 
  (float2) (1,0),
  (float2) (0.980785f,0.19509f), 
  (float2) (0.92388f,0.382683f),
  (float2) (0.83147f,0.55557f),
  (float2) (0.707107f,0.707107f),
  (float2) (0.55557f,0.83147f),
  (float2) (0.382683f,0.92388f),
  (float2) (0.19509f,0.980785f),
  (float2) (0,1),
  (float2) (-0.19509f,0.980785f),
  (float2) (-0.382683f,0.92388f),
  (float2) (-0.55557f,0.83147f),
  (float2) (-0.707107f,0.707107f),
  (float2) (-0.83147f,0.55557f),
  (float2) (-0.92388f,0.382683f),
  (float2) (-0.980785f,0.19509f),
};


inline uint5_t bit_reverse_32(uint5_t n)
{
  return ((n & 0x01) << 4) | ((n & 0x02) << 2) | ((n & 0x04)) | ((n & 0x08) >> 2) | ((n & 0x10) >> 4);
}


inline uint1_t parity(uint5_t n)
{
  uint1_t p = 0;

#pragma unroll
  for (unsigned i = 0; i < 5; i ++)
    p ^= (n >> i);

  return p;
}


inline void fft_32x32_4(/*float8 out[restrict 32][32],*/ float8 out_upper[restrict 16][32], float8 out_lower[restrict 16][32],  const float8 in[restrict 32][32], int sign)
{
 float8 left[32][16], right[32][16];

#pragma loop_coalesce 2
#pragma ii 1
  for (uint2_t dim = 0; dim < 2; dim ++) {
#pragma ii 1
    for (uint3_t level = 0; level < 5; level ++) {
#pragma loop_coalesce 2
#pragma ii 1
#pragma ivdep
      for (uint5_t k = 0; k < 16; k ++) {
        uint5_t mask = (1 << level) - 1;
        uint5_t e_index = bit_reverse_32(((k & ~mask) << 1) | (k & mask));
        uint5_t o_index = e_index | (16 >> level);

        uint4_t w_index = k << (4 - level);
        float2  w       = cosisin_16[w_index];

#pragma ivdep
        for (uint6_t n = 0; n < 32; n ++) {
          uint1_t swap    = parity(e_index ^ n);
          uint5_t l_index = swap ? o_index : e_index;
          uint5_t r_index = swap ? e_index : o_index;

          float8 l = dim == 0 ?  left[n][l_index / 2] :  left[l_index][n / 2];
          float8 r = dim == 0 ? right[n][r_index / 2] : right[r_index][n / 2];
          float8 e = swap ? r : l;
          float8 o = swap ? l : r;

          if (dim == 0 && level == 0) {
	    uint5_t e_index = bit_reverse_32(k << 1);
	    uint5_t o_index = e_index | 16;
	    e = in[n][e_index];
	    o = in[n][o_index];
          }

          float8 sum;
          sum.even = w.x * o.even - w.y * o.odd + e.even;
          sum.odd  = w.y * o.even + w.x * o.odd + e.odd;
          float8 diff = 2 * e - sum;

          * (dim == 0 ? & left[n][l_index / 2] : & left[l_index][n / 2]) = swap ? diff : sum;
          * (dim == 0 ? &right[n][r_index / 2] : &right[r_index][n / 2]) = swap ? sum  : diff;

          if (dim == 1 && level == 4) {
#if 0
	    uint5_t e_index = k;
	    uint5_t o_index = e_index | 16;
	    uint5_t n_index = bit_reverse_32(n);
	    out[e_index][n_index] = sum;
	    out[o_index][n_index] = diff;
#else
	    uint5_t n_index = bit_reverse_32(n);
	    out_lower[k][n_index] = sum;
	    out_upper[k][n_index] = diff;
#endif
	  }
        }
      }
    }
  }
}


void do_fft(__global float2 out[restrict NR_POLARIZATIONS][NR_PIXELS], float8 data[SUBGRID_SIZE][SUBGRID_SIZE])
{
#if 0
  float8 tmp[SUBGRID_SIZE][SUBGRID_SIZE] __attribute__((bank_bits(9))); // compiler does not understand how to do this without arbiter

  fft_32x32_4(tmp, data, FFT_BACKWARD);
#else
  float8 tmp_upper[SUBGRID_SIZE / 2][SUBGRID_SIZE], tmp_lower[SUBGRID_SIZE / 2][SUBGRID_SIZE];

  fft_32x32_4(tmp_upper, tmp_lower, data, FFT_BACKWARD);
#endif
  
#pragma loop_coalesce
#pragma ii 1
  for (unsigned short polarization = 0; polarization < NR_POLARIZATIONS; polarization ++)
    for (unsigned short pixel = 0; pixel < NR_PIXELS; pixel ++) {
      unsigned short x = pixel % SUBGRID_SIZE, y = pixel / SUBGRID_SIZE;
      float8         tmp = y >= 16 ? tmp_upper[y & 15][x] : tmp_lower[y][x];

      out[polarization][pixel] = (float2) (tmp[2 * polarization + REAL], tmp[2 * polarization + IMAG]);
    }
}


void post_process(
  __global float2 output_subgrid[restrict NR_POLARIZATIONS][NR_PIXELS],
  float8 subgrids[restrict NR_GRIDDERS][NR_PIXELS],
  __global const volatile float8 aterms[restrict NR_STATIONS][NR_PIXELS],
  uint2 stationPair,
  const float spheroidal[NR_PIXELS],
  unsigned short gridder)
{
  float8 subgrid2[NR_PIXELS];

  apply_aterm_and_spheroidal(subgrid2, subgrids, aterms, stationPair, spheroidal, gridder);
  do_fft(output_subgrid, (void *) subgrid2);
}


__attribute__((max_global_work_dim(0)))
__kernel void gridder(
  __global float2 output_subgrids[restrict][NR_POLARIZATIONS][NR_PIXELS],
  __global const volatile float8 visibilities[restrict][NR_TIMESTEPS][NR_CHANNELS],
  __global const volatile float wave_numbers[restrict NR_CHANNELS],
  __global const volatile float lmn[restrict NR_PIXELS][3],
  __global const volatile float uvw[restrict][NR_TIMESTEPS][3],
  __global const volatile float uvw_offsets[restrict][3],
  __global const volatile float8 aterms[restrict NR_STATIONS][NR_PIXELS],
  __global const volatile uint2 stationPairs[restrict],
  __global const volatile float spheroidal[restrict NR_PIXELS],
  unsigned nr_subgrids)
{
  float _wave_numbers[NR_CHANNELS];

#pragma ii 2
  for (unsigned short chan = 0; chan < NR_CHANNELS; chan ++)
    _wave_numbers[chan] = wave_numbers[chan] * cosisin_phase_correction_factor();

  float3 _lmn[NR_PIXELS];

#pragma ii 2
  for (unsigned short pixel = 0; pixel < NR_PIXELS; pixel ++)
    _lmn[pixel] = load_float3(lmn[pixel]);

  float _spheroidal[NR_PIXELS];

#pragma ii 2
  for (unsigned short pixel = 0; pixel < NR_PIXELS; pixel ++)
    _spheroidal[pixel] = spheroidal[pixel];

#pragma max_concurrency 8
#pragma ivdep array(output_subgrids)
#pragma ii 2
  for (unsigned subgrid_major = 0; subgrid_major < nr_subgrids; subgrid_major += NR_GRIDDERS) {
    visibilities[0][0][0]; // dummy statement that keeps the 20.3 compiler from crashing

    float8 subgrids[NR_GRIDDERS][NR_PIXELS];
    float8 pixels[NR_GRIDDERS][NR_PIXELS];
    float phase_offsets[NR_GRIDDERS][NR_PIXELS];

#if NR_GRIDDERS == 2 // work around compiler bug???
#pragma ii 2
    for (unsigned short gridder = 0; gridder < NR_GRIDDERS; gridder ++) {
      // at the end of the outer loop, reprocess subgrid zero redundantly
      // if (subgrid_major + gridder >= nr_subgrids), instead of conditionally calling
      // do_subgrid() and post_process();
      unsigned subgrid = subgrid_major + gridder < nr_subgrids ? subgrid_major + gridder : 0;

      float3 uvw_offset = load_float3(uvw_offsets[subgrid]) * cosisin_phase_correction_factor();

      for (unsigned short pixel = 0; pixel < NR_PIXELS; pixel ++)
	phase_offsets[gridder][pixel] = dot(_lmn[pixel], uvw_offset);
    }
#else
    float3 _uvw_offsets[NR_GRIDDERS];

#pragma ii 2
    for (unsigned short gridder = 0; gridder < NR_GRIDDERS; gridder ++) {
      unsigned subgrid = subgrid_major + gridder < nr_subgrids ? subgrid_major + gridder : 0;

      _uvw_offsets[gridder] = load_float3(uvw_offsets[subgrid]) * cosisin_phase_correction_factor();
    }

#pragma ii 2
    for (unsigned short pixel = 0; pixel < NR_PIXELS; pixel ++) {
      float _phase_offsets[NR_GRIDDERS];

#pragma ii 2
      for (unsigned short gridder = 0; gridder < NR_GRIDDERS; gridder ++)
	_phase_offsets[gridder] = dot(_lmn[pixel], _uvw_offsets[gridder]);

#pragma unroll
      for (unsigned short gridder = 0; gridder < NR_GRIDDERS; gridder ++)
	phase_offsets[gridder][pixel] = _phase_offsets[gridder];
    }
#endif

#pragma max_concurrency 8
#pragma ivdep array(pixels)
#pragma ii 2
    for (unsigned short time_major = 0; time_major < NR_TIMESTEPS; time_major += 16) {
      float8 _visibilities[NR_GRIDDERS][16][NR_CHANNELS];
      float3 _uvw[NR_GRIDDERS][16];

#pragma ii 2
      for (unsigned short gridder = 0; gridder < NR_GRIDDERS; gridder ++) {
	unsigned subgrid = subgrid_major + gridder < nr_subgrids ? subgrid_major + gridder : 0;

	for (unsigned short time_minor = 0; time_minor < 16; time_minor ++)
	  _uvw[gridder][time_minor] = load_float3(uvw[subgrid][time_major + time_minor]);
      }

#pragma ii 2
      for (unsigned short gridder = 0; gridder < NR_GRIDDERS; gridder ++) {
	unsigned subgrid = subgrid_major + gridder < nr_subgrids ? subgrid_major + gridder : 0;

#pragma loop_coalesce
	for (unsigned short time_minor = 0; time_minor < 16; time_minor ++)
	  for (unsigned short chan = 0; chan < NR_CHANNELS; chan ++)
	    _visibilities[gridder][time_minor][chan] = visibilities[subgrid][time_major + time_minor][chan];
      }

#pragma ivdep array(pixels)
#pragma ii 2
      for (unsigned short time_minor = 0; time_minor < 16; time_minor ++) {
#pragma unroll
	for (unsigned short gridder = 0; gridder < NR_GRIDDERS; gridder ++) {
	  unsigned subgrid = subgrid_major + gridder < nr_subgrids ? subgrid_major + gridder : 0;

#if 0
	  do_subgrid(subgrids[gridder], pixels[gridder], _visibilities[gridder][time_minor], _wave_numbers, _lmn, _uvw[gridder][time_minor], phase_offsets[gridder], time_major + time_minor);
#else
	  float3 tmp_uvws[NR_GRIDDERS] __attribute__((register));

#pragma unroll
	  for (unsigned const_gridder = 0; const_gridder < NR_GRIDDERS; const_gridder ++)
	    tmp_uvws[const_gridder] = _uvw[const_gridder][time_minor];

	  do_subgrid(subgrids[gridder], pixels[gridder], _visibilities[gridder][time_minor], _wave_numbers, _lmn, tmp_uvws[gridder], phase_offsets[gridder], time_major + time_minor);
#endif
	}
      }
    }

#pragma ii 2
    for (unsigned short gridder = 0; gridder < NR_GRIDDERS; gridder ++) {
      unsigned subgrid = subgrid_major + gridder < nr_subgrids ? subgrid_major + gridder : 0;

      post_process(output_subgrids[subgrid], subgrids, aterms, stationPairs[subgrid], _spheroidal, gridder);
    }
  }
}
