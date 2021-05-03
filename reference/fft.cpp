#include "fft.h"

#define NR_POLARIZATIONS 4

void kernel_fft(
	long size,
	long batch,
    fftwf_complex *data,
    int sign    // Note: -1=FFTW_FORWARD, 1=FFTW_BACKWARD
	)
{

    // 2D FFT
    int rank = 2;

    // For grids of size*size elements
    int n[] = {(int) size, (int) size};

    // Set stride
    int istride = 1;
    int ostride = istride;

    // Set dist
    int idist = n[0] * n[1];
    int odist = idist;

    // Planner flags
    int flags = FFTW_ESTIMATE;

    // Create plan
    fftwf_plan plan;
    #pragma omp critical
    {
        plan = fftwf_plan_many_dft(
        rank, n, batch * NR_POLARIZATIONS, data, n,
        istride, idist, data, n,
        ostride, odist, sign, flags);
    }

    // Execute FFTs
    fftwf_execute_dft(plan, data, data);

    // Destroy plan
    #pragma omp critical
    fftwf_destroy_plan(plan);
}
