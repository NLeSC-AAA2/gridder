#include <fftw3.h>

void kernel_fft(
	long size,
	long batch,
    fftwf_complex *data,
    int sign);
