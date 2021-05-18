#include <complex>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#define __CL_ENABLE_EXCEPTIONS
#include <CL/cl.hpp>


std::string errorMessage(cl_int error)
{
  switch (error) {
  case CL_SUCCESS:                            return "Success!";
  case CL_DEVICE_NOT_FOUND:                   return "Device not found.";
  case CL_DEVICE_NOT_AVAILABLE:               return "Device not available";
  case CL_COMPILER_NOT_AVAILABLE:             return "Compiler not available";
  case CL_MEM_OBJECT_ALLOCATION_FAILURE:      return "Memory object allocation failure";
  case CL_OUT_OF_RESOURCES:                   return "Out of resources";
  case CL_OUT_OF_HOST_MEMORY:                 return "Out of host memory";
  case CL_PROFILING_INFO_NOT_AVAILABLE:       return "Profiling information not available";
  case CL_MEM_COPY_OVERLAP:                   return "Memory copy overlap";
  case CL_IMAGE_FORMAT_MISMATCH:              return "Image format mismatch";
  case CL_IMAGE_FORMAT_NOT_SUPPORTED:         return "Image format not supported";
  case CL_BUILD_PROGRAM_FAILURE:              return "Program build failure";
  case CL_MAP_FAILURE:                        return "Map failure";
  case CL_INVALID_VALUE:                      return "Invalid value";
  case CL_INVALID_DEVICE_TYPE:                return "Invalid device type";
  case CL_INVALID_PLATFORM:                   return "Invalid platform";
  case CL_INVALID_DEVICE:                     return "Invalid device";
  case CL_INVALID_CONTEXT:                    return "Invalid context";
  case CL_INVALID_QUEUE_PROPERTIES:           return "Invalid queue properties";
  case CL_INVALID_COMMAND_QUEUE:              return "Invalid command queue";
  case CL_INVALID_HOST_PTR:                   return "Invalid host pointer";
  case CL_INVALID_MEM_OBJECT:                 return "Invalid memory object";
  case CL_INVALID_IMAGE_FORMAT_DESCRIPTOR:    return "Invalid image format descriptor";
  case CL_INVALID_IMAGE_SIZE:                 return "Invalid image size";
  case CL_INVALID_SAMPLER:                    return "Invalid sampler";
  case CL_INVALID_BINARY:                     return "Invalid binary";
  case CL_INVALID_BUILD_OPTIONS:              return "Invalid build options";
  case CL_INVALID_PROGRAM:                    return "Invalid program";
  case CL_INVALID_PROGRAM_EXECUTABLE:         return "Invalid program executable";
  case CL_INVALID_KERNEL_NAME:                return "Invalid kernel name";
  case CL_INVALID_KERNEL_DEFINITION:          return "Invalid kernel definition";
  case CL_INVALID_KERNEL:                     return "Invalid kernel";
  case CL_INVALID_ARG_INDEX:                  return "Invalid argument index";
  case CL_INVALID_ARG_VALUE:                  return "Invalid argument value";
  case CL_INVALID_ARG_SIZE:                   return "Invalid argument size";
  case CL_INVALID_KERNEL_ARGS:                return "Invalid kernel arguments";
  case CL_INVALID_WORK_DIMENSION:             return "Invalid work dimension";
  case CL_INVALID_WORK_GROUP_SIZE:            return "Invalid work group size";
  case CL_INVALID_WORK_ITEM_SIZE:             return "Invalid work item size";
  case CL_INVALID_GLOBAL_OFFSET:              return "Invalid global offset";
  case CL_INVALID_EVENT_WAIT_LIST:            return "Invalid event wait list";
  case CL_INVALID_EVENT:                      return "Invalid event";
  case CL_INVALID_OPERATION:                  return "Invalid operation";
  case CL_INVALID_GL_OBJECT:                  return "Invalid OpenGL object";
  case CL_INVALID_BUFFER_SIZE:                return "Invalid buffer size";
  case CL_INVALID_MIP_LEVEL:                  return "Invalid mip-map level";
  case CL_INVALID_GLOBAL_WORK_SIZE:           return "Invalid global work size";
#if defined CL_INVALID_PROPERTY
  case CL_INVALID_PROPERTY:                   return "Invalid property";
#endif
#if defined CL_INVALID_IMAGE_DESCRIPTOR
  case CL_INVALID_IMAGE_DESCRIPTOR:           return "Invalid image descriptor";
#endif
#if defined CL_INVALID_COMPILER_OPTIONS
  case CL_INVALID_COMPILER_OPTIONS:           return "Invalid compiler options";
#endif
#if defined CL_INVALID_LINKER_OPTIONS
  case CL_INVALID_LINKER_OPTIONS:             return "Invalid linker options";
#endif
#if defined CL_INVALID_DEVICE_PARTITION_COUNT
  case CL_INVALID_DEVICE_PARTITION_COUNT:     return "Invalid device partition count";
#endif
  default:                                    std::stringstream str;
    str << "Unknown (" << error << ')';
    return str.str();
  }
}


void createContext(cl::Context &context, std::vector<cl::Device> &devices)
{
  const char *platformName = getenv("PLATFORM");

#if defined __linux__
  if (platformName == 0)
#endif
  platformName = "Intel(R) FPGA SDK for OpenCL(TM)";
  //platformName = Intel(R) OpenCL";
  //platformName = "NVIDIA CUDA";
  //platformName = "AMD Accelerated Parallel Processing";

  cl_device_type type = CL_DEVICE_TYPE_DEFAULT;

  const char *deviceType = getenv("TYPE");

  if (deviceType != 0) {
    if (strcmp(deviceType, "GPU") == 0)
      type = CL_DEVICE_TYPE_GPU;
    else if (strcmp(deviceType, "CPU") == 0)
      type = CL_DEVICE_TYPE_CPU;
    else if (strcmp(deviceType, "ACCELERATOR") == 0)
      type = CL_DEVICE_TYPE_ACCELERATOR;
    else
      std::cerr << "Unrecognized device type: " << deviceType;
  }

  const char *deviceName = getenv("DEVICE");

  std::vector<cl::Platform> platforms;
  cl::Platform::get(&platforms);

  for (cl::Platform &platform : platforms) {
    std::clog << "Platform name: " << platform.getInfo<CL_PLATFORM_NAME>() << std::endl;
    std::clog << "Platform version: " << platform.getInfo<CL_PLATFORM_VERSION>() << std::endl;
    std::clog << "Platform extensions: " << platform.getInfo<CL_PLATFORM_EXTENSIONS>() << std::endl;
  }

  for (cl::Platform &platform : platforms) {
    if (strcmp(platform.getInfo<CL_PLATFORM_NAME>().c_str(), platformName) == 0) {
      platform.getDevices(type, &devices);

      for (cl::Device &device : devices) {
	std::clog << "Device: " << device.getInfo<CL_DEVICE_NAME>() << ", "
		      "mem: " << device.getInfo<CL_DEVICE_GLOBAL_MEM_SIZE>() / (1024 * 1024) << " MB, max alloc: " << device.getInfo<CL_DEVICE_MAX_MEM_ALLOC_SIZE>() / (1024 * 1024) << " MB"
//#if defined CL_DEVICE_TOPOLOGY_AMD
//		      ", PCIe bus 0x" << std::hex << std::setiosflags (std::ios::uppercase) << std::setfill('0') << std::setw(2) << (device.getInfo<CL_DEVICE_TOPOLOGY_AMD>().pcie.bus & 255) << std::dec
//#endif
	          << std::endl;
      }

      context = cl::Context(devices);
      return;
    }
  }

  std::cerr << "Platform not found: \"" << platformName << '"' << std::endl;
  exit(1);
}


cl::Program createProgramFromSources(cl::Context &context, std::vector<cl::Device> &devices, const std::string &sources, const char *args)
{
  std::stringstream cmd;
  cmd << "#include \"" << sources << '"' << std::endl;
#if defined CL_VERSION_1_2
  cl::Program program(context, cmd.str());
#else
  std::string str = cmd.str();
  cl::Program::Sources src(1, std::make_pair(str.c_str(), str.length()));
  cl::Program program(context, src);
#endif

  try {
    program.build(devices, args);
    std::string msg;
    program.getBuildInfo(devices[0], CL_PROGRAM_BUILD_LOG, &msg);

    std::clog << msg << std::endl;
  } catch (cl::Error &error) {
    if (strcmp(error.what(), "clBuildProgram") == 0) {
      std::string msg;
      program.getBuildInfo(devices[0], CL_PROGRAM_BUILD_LOG, &msg);

      std::cerr << msg << std::endl;
      exit(1);
    } else {
      throw;
    }
  }

#if 0
  std::vector<size_t> binarySizes = program.getInfo<CL_PROGRAM_BINARY_SIZES>();
  std::vector<char *> binaries = program.getInfo<CL_PROGRAM_BINARIES>();

  for (unsigned i = 0; i < binaries.size(); i++) {
    std::stringstream filename;
    filename << sources << '-' << i << ".ptx";
    std::ofstream(filename.str().c_str(), std::ofstream::binary).write(binaries[i], binarySizes[i]);
  }

  for (unsigned b = 0; b < binaries.size(); b++)
    delete [] binaries[b];
#endif

  return program;
}


cl::Program createProgramFromBinaries(cl::Context &context, std::vector<cl::Device> &devices, const std::string &name) 
{
  std::ifstream ifs(name, std::ios::in | std::ios::binary);
  std::string str((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());

  cl::Program::Binaries binaries(devices.size(), std::make_pair(str.c_str(), str.length()));
  return cl::Program(context, devices, binaries);
} 


#define NR_RECEIVERS            576
#define NR_SAMPLES_PER_CHANNEL  1024
#define NR_CHANNELS             64
#define NR_TAPS                 16

typedef signed char int8_t;

typedef std::complex<int8_t> InputType[NR_RECEIVERS][NR_SAMPLES_PER_CHANNEL + NR_TAPS - 1][NR_CHANNELS];
typedef float FIR_FilterWeightsType[NR_CHANNELS][NR_TAPS];
typedef float BandPassWeightsType[NR_CHANNELS];
typedef struct { float atBegin, afterEnd; } DelaysType[NR_RECEIVERS];
typedef std::complex<float> OutputType[NR_RECEIVERS][NR_SAMPLES_PER_CHANNEL][NR_CHANNELS];

#if !defined CL_VERSION_1_2
#define CL_MEM_HOST_READ_ONLY	0
#define CL_MEM_HOST_WRITE_ONLY	0
#endif



int main(int argc, char * argv[])
{
  if ( argc != 2) {
    std::cerr << "Usage: " << argv[0] << " <AOCX file>" << std::endl;
    exit(-1);
  }
  try {
    cl::Context             context;
    std::vector<cl::Device> devices;
    createContext(context, devices);

    cl::Program program(createProgramFromBinaries(context, devices, std::string(argv[1])));
    cl::CommandQueue queue(context, devices[0], CL_QUEUE_PROFILING_ENABLE);

    cl::Buffer inputBuffer(context, CL_MEM_READ_ONLY | CL_MEM_HOST_WRITE_ONLY | CL_MEM_ALLOC_HOST_PTR, sizeof(InputType));
    InputType &input = * (InputType *) queue.enqueueMapBuffer(inputBuffer, CL_TRUE, CL_MAP_WRITE, 0, sizeof(InputType));
    memset(input, 0, sizeof(InputType));
    input[4][94][5] = std::complex<float>(4,0);
    input[4][95][5] = std::complex<float>(5,0);
    queue.enqueueUnmapMemObject(inputBuffer, input);

    cl::Buffer firFilterWeightsBuffer(context, CL_MEM_READ_ONLY | CL_MEM_HOST_WRITE_ONLY | CL_MEM_ALLOC_HOST_PTR, sizeof(FIR_FilterWeightsType));
    FIR_FilterWeightsType &firFilterWeights = * (FIR_FilterWeightsType *) queue.enqueueMapBuffer(firFilterWeightsBuffer, CL_TRUE, CL_MAP_WRITE, 0, sizeof(FIR_FilterWeightsType));
    memset(firFilterWeights, 0, sizeof(FIR_FilterWeightsType));
    firFilterWeights[5][3] = 2;
    firFilterWeights[5][4] = 3;
    queue.enqueueUnmapMemObject(firFilterWeightsBuffer, firFilterWeights);

    cl::Buffer bandPassWeightsBuffer(context, CL_MEM_READ_ONLY | CL_MEM_HOST_WRITE_ONLY | CL_MEM_ALLOC_HOST_PTR, sizeof(BandPassWeightsType));
    BandPassWeightsType &bandPassWeights = * (BandPassWeightsType *) queue.enqueueMapBuffer(bandPassWeightsBuffer, CL_TRUE, CL_MAP_WRITE, 0, sizeof(BandPassWeightsType));

    for (int ch = 0; ch < NR_CHANNELS; ch ++)
      bandPassWeights[ch] = ch == 5 ? 2 : 9;

    queue.enqueueUnmapMemObject(bandPassWeightsBuffer, bandPassWeights);

    cl::Buffer delaysBuffer(context, CL_MEM_READ_ONLY | CL_MEM_HOST_WRITE_ONLY | CL_MEM_ALLOC_HOST_PTR, sizeof(DelaysType));
    DelaysType &delays = * (DelaysType *) queue.enqueueMapBuffer(delaysBuffer, CL_TRUE, CL_MAP_WRITE, 0, sizeof(DelaysType));
    memset(delays, 0, sizeof(DelaysType));
    delays[4].atBegin = 1e-6;
    delays[4].afterEnd = 1.1e-6;
    queue.enqueueUnmapMemObject(delaysBuffer, delays);

    cl::Buffer outputBuffer(context, CL_MEM_WRITE_ONLY | CL_MEM_HOST_READ_ONLY | CL_MEM_ALLOC_HOST_PTR, sizeof(OutputType));

    cl::Event event;

#if 1
    cl::Kernel kernel(program, "fused");
    kernel.setArg(0, inputBuffer);
    kernel.setArg(1, firFilterWeightsBuffer);
    kernel.setArg(2, 60e6f);
    kernel.setArg(3, delaysBuffer);
    kernel.setArg(4, bandPassWeightsBuffer);
    kernel.setArg(5, outputBuffer);

    queue.enqueueTask(kernel, 0, &event);
#else
    cl::Kernel readerKernel(program, "reader");
    readerKernel.setArg(0, inputBuffer);
    cl::CommandQueue(context, devices[0], CL_QUEUE_PROFILING_ENABLE).enqueueTask(readerKernel);

    cl::Kernel firFilterKernel(program, "firFilterKernel");
    firFilterKernel.setArg(0, firFilterWeightsBuffer);
    cl::CommandQueue(context, devices[0], CL_QUEUE_PROFILING_ENABLE).enqueueTask(firFilterKernel);

    cl::Kernel fftKernel(program, "fftKernel");
    cl::CommandQueue(context, devices[0], CL_QUEUE_PROFILING_ENABLE).enqueueTask(fftKernel);

    cl::Kernel delayAndBandPassKernel(program, "delayAndBandPassKernel");
    delayAndBandPassKernel.setArg(0, 60e6f);
    delayAndBandPassKernel.setArg(1, delaysBuffer);
    delayAndBandPassKernel.setArg(2, bandPassWeightsBuffer);
    delayAndBandPassKernel.setArg(3, outputBuffer);
    cl::CommandQueue delayAndBandPassQueue(context, devices[0], CL_QUEUE_PROFILING_ENABLE);
    delayAndBandPassQueue.enqueueTask(delayAndBandPassKernel, 0, &event);
    delayAndBandPassQueue.finish();
#endif

#if 1
    event.wait();
    cl_ulong start = event.getProfilingInfo<CL_PROFILING_COMMAND_START>();
    cl_ulong stop  = event.getProfilingInfo<CL_PROFILING_COMMAND_END>();
    double seconds = (stop - start) / 1e9;
    std::cout << "runtime = " << seconds << " s" << std::endl;
#endif

    OutputType &output = * (OutputType *) queue.enqueueMapBuffer(outputBuffer, CL_TRUE, CL_MAP_READ, 0, sizeof(OutputType));

    for (int recv = 0; recv < NR_RECEIVERS; recv ++)
      for (int time = 0; time < NR_SAMPLES_PER_CHANNEL; time ++)
	for (int ch = 0; ch < NR_CHANNELS; ch ++)
	  if (output[recv][time][ch] != std::complex<float>(0))
	      std::cout << "output[" << recv << "][" << time << "][" << ch << "] = " << output[recv][time][ch] << std::endl;

    queue.enqueueUnmapMemObject(outputBuffer, output);
   } catch (cl::Error &error) {
     std::cerr << "caught cl::Error: " << error.what() << ": " << errorMessage(error.err()) << std::endl;
   }

  return 0;
}
