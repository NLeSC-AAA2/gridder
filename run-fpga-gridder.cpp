#include <fftw3.h>
#include <omp.h>

#include <CL/cl_ext_intelfpga.h>

#include "common/common.h"
#include "common/print.h"
#include "common/init.h"
#include "reference/fft.h"
#include "reference/gridder.h"

int main(int argc, char **argv)
{
    if (argc > 2) {
        std::cerr << "usage: " << argv[0] << " [gridder.aocx]" << std::endl;
        exit(1);
    }

    // Initialize OpenCL
    cl::Context context;
    std::vector<cl::Device> devices;
    init(context, devices);
    cl::Device &device = devices[0];

    // Get program
    std::string filename_bin = std::string(argc == 2 ? argv[1] : "gridder.aocx");
    cl::Program program = get_program(context, device, filename_bin);

    // Allocate data structures on host
    std::clog << ">>> Initialize data structures on host" << std::endl;
    idg::Array2D<idg::UVWCoordinate<float>> uvw(NR_BASELINES, NR_TIMESTEPS);
    idg::Array3D<idg::Visibility<std::complex<float>>> visibilities(NR_BASELINES, NR_TIMESTEPS, NR_CHANNELS);
    idg::Array1D<idg::Baseline> baselines(NR_BASELINES);
    idg::Array4D<idg::Matrix2x2<std::complex<float>>> aterms(NR_TIMESLOTS, NR_STATIONS, SUBGRID_SIZE, SUBGRID_SIZE);
    idg::Array1D<float> frequencies(NR_CHANNELS);
    idg::Array1D<float> wavenumbers(NR_CHANNELS);
    idg::Array2D<float> spheroidal(SUBGRID_SIZE, SUBGRID_SIZE);
    idg::Array4D<std::complex<float>> subgrids(NR_SUBGRIDS, NR_CORRELATIONS, SUBGRID_SIZE, SUBGRID_SIZE);
    idg::Array1D<idg::Metadata> metadata(NR_SUBGRIDS);

    // Initialize random number generator
    srand(0);

    // Initialize data structures
    initialize_uvw(GRID_SIZE, uvw);
    initialize_frequencies(frequencies);
    initialize_wavenumbers(frequencies, wavenumbers);
    initialize_visibilities(GRID_SIZE, IMAGE_SIZE, frequencies, uvw, visibilities);
    initialize_baselines(NR_STATIONS, baselines);
    initialize_spheroidal(spheroidal);
    initialize_aterms(spheroidal, aterms);
    initialize_metadata(GRID_SIZE, NR_TIMESLOTS, NR_TIMESTEPS_SUBGRID, baselines, metadata);

    // Run reference
    std::clog << ">>> Run reference" << std::endl;
    idg::Array4D<std::complex<float>> subgrids_ref(NR_SUBGRIDS, NR_CORRELATIONS, SUBGRID_SIZE, SUBGRID_SIZE);
#if 1 
    kernel_gridder(
        NR_SUBGRIDS, GRID_SIZE, SUBGRID_SIZE, IMAGE_SIZE, W_STEP, NR_CHANNELS, NR_STATIONS,
        uvw.data(), wavenumbers.data(), (std::complex<float> *) visibilities.data(),
        (float *) spheroidal.data(), (std::complex<float> *) aterms.data(), metadata.data(),
        subgrids_ref.data());
    kernel_fft(SUBGRID_SIZE, NR_SUBGRIDS, (fftwf_complex *) subgrids_ref.data(), 1);

    // Print subgrids
    //print_subgrid(subgrids_ref, 0);
    print_subgrid(subgrids_ref, NR_SUBGRIDS - 1);
#endif

    // Additional data structures for FPGA
    idg::Array3D<float> lmn_fpga(SUBGRID_SIZE, SUBGRID_SIZE, 3);
    idg::Array2D<float> uvw_offsets_fpga(NR_SUBGRIDS, 3);
    idg::Array3D<idg::Visibility<std::complex<float>>> visibilities_fpga(NR_SUBGRIDS, NR_TIMESTEPS_SUBGRID, NR_CHANNELS);
    idg::Array2D<idg::UVWCoordinate<float>> uvw_fpga(NR_SUBGRIDS, NR_TIMESTEPS_SUBGRID);
    idg::Array1D<idg::Baseline> baselines_fpga(NR_SUBGRIDS);

    // Initialize data structures
    initialize_lmn(IMAGE_SIZE, lmn_fpga);
    initialize_uvw_offsets(SUBGRID_SIZE, GRID_SIZE, IMAGE_SIZE, W_STEP, metadata, uvw_offsets_fpga);

    // Transpose data structures
    for (unsigned bl = 0; bl < NR_BASELINES; bl++) {
        for (unsigned ts = 0; ts < NR_TIMESLOTS; ts++) {
            unsigned s = bl * NR_TIMESLOTS + ts;

            baselines_fpga(s) = baselines(bl);

            for (unsigned t = 0; t < NR_TIMESTEPS_SUBGRID; t++) {
                unsigned time = ts * NR_TIMESTEPS_SUBGRID + t;

                uvw_fpga(s, t) = uvw(bl, time);

                for (unsigned chan = 0; chan < NR_CHANNELS; chan++) {
                    visibilities_fpga(s, t, chan) = visibilities(bl, time, chan);
                }
            }
        }
    }

    // Setup command queues
    std::vector<cl::CommandQueue> queues(8);

    for (cl::CommandQueue &queue : queues) {
        queue = cl::CommandQueue(context, device, CL_QUEUE_PROFILING_ENABLE);
    }

    // Allocate data structures on device
    std::clog << ">>> Allocate data structures on device" << std::endl;
    cl::Buffer d_uvw          = cl::Buffer(context, CL_MEM_READ_ONLY, uvw_fpga.bytes());
    cl::Buffer d_wavenumbers  = cl::Buffer(context, CL_MEM_READ_ONLY, wavenumbers.bytes());
    cl::Buffer d_visibilities = cl::Buffer(context, CL_MEM_READ_ONLY, visibilities_fpga.bytes());
    cl::Buffer d_spheroidal   = cl::Buffer(context, CL_MEM_READ_ONLY, spheroidal.bytes());
    cl::Buffer d_aterms       = cl::Buffer(context, CL_MEM_READ_ONLY, aterms.bytes());
    cl::Buffer d_subgrids     = cl::Buffer(context, CL_MEM_READ_WRITE, subgrids.bytes());
    cl::Buffer d_lmn          = cl::Buffer(context, CL_MEM_READ_ONLY, lmn_fpga.bytes());
    cl::Buffer d_uvw_offsets  = cl::Buffer(context, CL_MEM_READ_ONLY, uvw_offsets_fpga.bytes());
    cl::Buffer d_baselines    = cl::Buffer(context, CL_MEM_READ_ONLY, baselines_fpga.bytes());
    cl::Buffer d_nr_subgrids;

    cl::Event computeDone;

#if 0
    std::vector<cl::Event> lmnCopied(1);

    // Copy data structures to device
    std::clog << ">>> Copy data structures to device" << std::endl;
    queues[0].enqueueWriteBuffer(d_uvw,          CL_TRUE, 0, uvw_fpga.bytes(), uvw_fpga.data());
    queues[1].enqueueWriteBuffer(d_lmn,          CL_TRUE, 0, lmn_fpga.bytes(), lmn_fpga.data(), nullptr, &lmnCopied[0]);
    queues[2].enqueueWriteBuffer(d_uvw_offsets,  CL_TRUE, 0, uvw_offsets_fpga.bytes(), uvw_offsets_fpga.data());
    queues[3].enqueueWriteBuffer(d_wavenumbers,  CL_TRUE, 0, wavenumbers.bytes(), wavenumbers.data());
    queues[4].enqueueWriteBuffer(d_visibilities, CL_TRUE, 0, visibilities_fpga.bytes(), visibilities_fpga.data());
    queues[5].enqueueWriteBuffer(d_aterms,       CL_TRUE, 0, aterms.bytes(), aterms.data());
    queues[5].enqueueWriteBuffer(d_baselines,    CL_TRUE, 0, baselines_fpga.bytes(), baselines_fpga.data());
    queues[6].enqueueWriteBuffer(d_spheroidal,   CL_TRUE, 0, spheroidal.bytes(), spheroidal.data());

    // Work around emulator bug
    if (getenv("CL_CONTEXT_EMULATOR_DEVICE_INTELFPGA")) {
        d_nr_subgrids = cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(unsigned));
        unsigned nr_subgrids = NR_SUBGRIDS;
        queues[0].enqueueWriteBuffer(d_nr_subgrids, CL_TRUE, 0, sizeof(unsigned), &nr_subgrids);
    }

    // Setup FPGA kernels
    cl::Kernel uvwReaderKernel(program, "uvw_reader");
    uvwReaderKernel.setArg(0, d_uvw);
    uvwReaderKernel.setArg(1, NR_SUBGRIDS);

    cl::Kernel lmnReaderKernel(program, "lmn_reader");
    lmnReaderKernel.setArg(0, d_lmn);
    lmnReaderKernel.setArg(1, NR_SUBGRIDS);

    cl::Kernel uvwOffsetsReaderKernel(program, "uvw_offsets_reader");
    uvwOffsetsReaderKernel.setArg(0, d_uvw_offsets);
    uvwOffsetsReaderKernel.setArg(1, d_lmn);
    uvwOffsetsReaderKernel.setArg(2, NR_SUBGRIDS);

    cl::Kernel waveNumbersReaderKernel(program, "wave_numbers_reader");
    waveNumbersReaderKernel.setArg(0, d_wavenumbers);

    if (getenv("CL_CONTEXT_EMULATOR_DEVICE_INTELFPGA")) {
        waveNumbersReaderKernel.setArg(1, d_nr_subgrids);
    } else {
        waveNumbersReaderKernel.setArg(1, NR_SUBGRIDS);
    }

    cl::Kernel visibilitiesReaderKernel(program, "visibilities_reader");
    visibilitiesReaderKernel.setArg(0, d_visibilities);
    visibilitiesReaderKernel.setArg(1, NR_SUBGRIDS);

    cl::Kernel aTermKernel(program, "aterm");
    aTermKernel.setArg(0, d_aterms);
    aTermKernel.setArg(1, d_baselines);
    aTermKernel.setArg(2, NR_SUBGRIDS);

    cl::Kernel spheroidalKernel(program, "spheroidal");
    spheroidalKernel.setArg(0, d_spheroidal);
    spheroidalKernel.setArg(1, NR_SUBGRIDS);

    cl::Kernel subgridWriterKernel(program, "subgrid_writer");
    subgridWriterKernel.setArg(0, d_subgrids);
    subgridWriterKernel.setArg(1, NR_SUBGRIDS);

    // Run FPGA kernels
    std::clog << ">>> Run fpga" << std::endl;
    try {
        queues[0].enqueueTask(uvwReaderKernel);
        queues[1].enqueueTask(lmnReaderKernel);
        queues[2].enqueueTask(uvwOffsetsReaderKernel, &lmnCopied);
        queues[3].enqueueTask(waveNumbersReaderKernel);
        queues[4].enqueueTask(visibilitiesReaderKernel);
        queues[5].enqueueTask(aTermKernel);
        queues[6].enqueueTask(spheroidalKernel);
        queues[7].enqueueTask(subgridWriterKernel, nullptr, &computeDone);
    } catch (cl::Error &error) {
        std::cerr << "Error launching kernel: " << error.what() << std::endl;
        exit(EXIT_FAILURE);
    }

    // Copy subgrid to host
    queues[7].enqueueReadBuffer(d_subgrids, CL_TRUE, 0, subgrids.bytes(), subgrids.data());
#else
    // Copy data structures to device
    std::clog << ">>> Copy data structures to device" << std::endl;
    queues[0].enqueueWriteBuffer(d_uvw,          CL_FALSE, 0, uvw_fpga.bytes(), uvw_fpga.data());
    queues[0].enqueueWriteBuffer(d_lmn,          CL_FALSE, 0, lmn_fpga.bytes(), lmn_fpga.data());
    queues[0].enqueueWriteBuffer(d_uvw_offsets,  CL_FALSE, 0, uvw_offsets_fpga.bytes(), uvw_offsets_fpga.data());
    queues[0].enqueueWriteBuffer(d_wavenumbers,  CL_FALSE, 0, wavenumbers.bytes(), wavenumbers.data());
    queues[0].enqueueWriteBuffer(d_visibilities, CL_FALSE, 0, visibilities_fpga.bytes(), visibilities_fpga.data());
    queues[0].enqueueWriteBuffer(d_aterms,       CL_FALSE, 0, aterms.bytes(), aterms.data());
    queues[0].enqueueWriteBuffer(d_baselines,    CL_FALSE, 0, baselines_fpga.bytes(), baselines_fpga.data());
    queues[0].enqueueWriteBuffer(d_spheroidal,   CL_FALSE, 0, spheroidal.bytes(), spheroidal.data());

    cl::Kernel gridderKernel(program, "gridder");
    gridderKernel.setArg(0, d_subgrids);
    gridderKernel.setArg(1, d_visibilities);
    gridderKernel.setArg(2, d_wavenumbers);
    gridderKernel.setArg(3, d_lmn);
    gridderKernel.setArg(4, d_uvw);
    gridderKernel.setArg(5, d_uvw_offsets);
    gridderKernel.setArg(6, d_aterms);
    gridderKernel.setArg(7, d_baselines);
    gridderKernel.setArg(8, d_spheroidal);
    gridderKernel.setArg(9, NR_SUBGRIDS);

    std::clog << ">>> Run fpga" << std::endl;

    try {
	queues[0].enqueueTask(gridderKernel, nullptr, &computeDone);
    } catch (cl::Error &error) {
	std::cerr << "Error launching kernel: " << error.what() << std::endl;
	exit(EXIT_FAILURE);
    }

    queues[0].enqueueReadBuffer(d_subgrids, CL_TRUE, 0, subgrids.bytes(), subgrids.data());
#endif

    // Print subgrids
    //print_subgrid(subgrids, 0);
    print_subgrid(subgrids, NR_SUBGRIDS - 1);

    // Compare subgrids
    std::clog << ">>> Difference between subgrids" << std::endl;
    //print_subgrid_diff(subgrids_ref, subgrids, 0);
    print_subgrid_diff(subgrids_ref, subgrids, NR_SUBGRIDS - 1);

    computeDone.wait();
    cl_ulong start = computeDone.getProfilingInfo<CL_PROFILING_COMMAND_START>();
    cl_ulong stop  = computeDone.getProfilingInfo<CL_PROFILING_COMMAND_END>();

    double milliseconds = (stop - start) / 1e6;
    std::cout << "runtime = " << milliseconds << " ms, " << std::endl;

#if 1
    {
	cl_int error;

	if (clGetProfileDataDeviceIntelFPGA(device(), program(), true, true, true, 0, nullptr, nullptr, &error) != CL_SUCCESS)
	    /*throw cl::Error(error, "clGetProfileDataDeviceIntelFPGA")*/; // ignore
    }
#endif

    return EXIT_SUCCESS;
}
