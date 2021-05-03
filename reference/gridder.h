#include <complex>

#include "types.h"

void kernel_gridder(
    const int   nr_subgrids,
    const int   gridsize,
    const int   subgridsize,
    const float imagesize,
    const float w_step_in_lambda,
    const int   nr_channels,
    const int   nr_stations,
    const idg::UVWCoordinate<float>* uvw,
    const float*                     wavenumbers,
    const std::complex<float>*       visibilities,
    const float*                     spheroidal,
    const std::complex<float>*       aterms,
    const idg::Metadata*             metadata,
          std::complex<float>*       subgrid);
