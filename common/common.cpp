#include "common.h"

using namespace std;

ostream &os = clog;

void init(
    cl::Context &context,
    vector<cl::Device> &devices)
{
    vector<cl::Platform> platforms;
    cl::Platform::get(&platforms);

    // The selected device
    int i = 0;
    const char *platform_name = getenv("PLATFORM");

    if (platform_name == 0)
      platform_name = getenv("CL_CONTEXT_EMULATOR_DEVICE_INTELFPGA") ? "Intel(R) FPGA Emulation Platform for OpenCL(TM)" : "Intel(R) FPGA SDK for OpenCL(TM)";

    os << ">>> OpenCL environment: " << endl;

    // Iterate all platforms
    for (cl::Platform &platform : platforms) {
        print_platform(platform);
	bool selected = platform.getInfo<CL_PLATFORM_NAME>() == platform_name;

        // Get devices for the current platform
        vector<cl::Device> devices_;
        platform.getDevices(CL_DEVICE_TYPE_ALL, &devices_);

        // Iterate all devices
        for (cl::Device &device : devices_) {
            if (selected)
                devices.push_back(device);

            print_device(device, selected);
            i++;
        }
    }
    os << endl;

    if (devices.size() == 0) {
        cerr << "Could not find any device in platform "  << platform_name << endl;
        exit(EXIT_FAILURE);
    }

    context = cl::Context(devices);
}

void print_platform(
    cl::Platform &platform)
{
    os << ">>> Platform: " << endl;
    os << "Name       : " << platform.getInfo<CL_PLATFORM_NAME>() << endl;
    os << "Version    : " << platform.getInfo<CL_PLATFORM_VERSION>() << endl;
    os << "Extensions : " << platform.getInfo<CL_PLATFORM_EXTENSIONS>() << endl;
    os <<  endl;
}

void print_device(
    cl::Device &device,
    bool marker)
{
    os << ">>> Device: ";
    if (marker) os << " (selected)";
    os << endl;
    os << "Name            : " << device.getInfo<CL_DEVICE_NAME>() << endl;
    os << "Driver version  : " << device.getInfo<CL_DRIVER_VERSION>() << endl;
    os << "Device version  : " << device.getInfo<CL_DEVICE_VERSION>() << endl;
    os << "Compute units   : " << device.getInfo<CL_DEVICE_MAX_COMPUTE_UNITS>() << endl;
    os << "Clock frequency : " << device.getInfo<CL_DEVICE_MAX_CLOCK_FREQUENCY>() << " MHz" << endl;
    os << "Global memory   : " << device.getInfo<CL_DEVICE_GLOBAL_MEM_SIZE>() * 1e-9 << " Gb" << endl;
    os << "Local memory    : " << device.getInfo<CL_DEVICE_LOCAL_MEM_SIZE>() * 1e-6 << " Mb" << endl;
    os << endl;
}

string get_source(
    string& filename)
{
    // Source directory
    string srcdir = "./cl";

    // All helper files to include in build
    vector<string> helper_files;
    helper_files.push_back("types.cl");
    helper_files.push_back("math.cl");

    // Store helper files in string
    stringstream source_helper_;

    for (int i = 0; i < helper_files.size(); i++) {
        // Get source filename
        stringstream source_file_name_;
        source_file_name_ << srcdir << "/" << helper_files[i];
        string source_file_name = source_file_name_.str();

        // Read source from file
        ifstream source_file(source_file_name.c_str());
        string source(istreambuf_iterator<char>(source_file),
                          (istreambuf_iterator<char>()));
        source_file.close();

        // Update source helper stream
        source_helper_ << source;
    }

    string source_helper = source_helper_.str();

    // Get source filename
    stringstream source_file_name_;
    source_file_name_ << srcdir << "/" << filename;
    string source_file_name = source_file_name_.str();

    // Read kernel source from file
    ifstream source_file(source_file_name.c_str());
    string source_kernel(
        istreambuf_iterator<char>(source_file),
        (istreambuf_iterator<char>()));
    source_file.close();

    // Construct full source file
    stringstream full_source;
    full_source << source_helper;
    full_source << source_kernel;

    return full_source.str();
}

string get_flags()
{
    return string("-cl-fast-relaxed-math");
}

void write_source(
    string& source,
    string& filename)
{
    cout << ">>> Writing source to: " << filename << endl
              << endl;
    ofstream source_output;
    source_output.open(filename, ofstream::out);
    source_output << source;
    source_output.close();
}

cl::Program get_program(
    cl::Context& context,
    cl::Device& device,
    string& filename)
{
    os << ">>> Loading program from binary: " << filename << endl;
    try {
        ifstream ifs(filename, ios::in | ios::binary);
        string str((istreambuf_iterator<char>(ifs)), istreambuf_iterator<char>());
        cl::Program::Binaries binaries(1, std::make_pair(str.c_str(), str.length()));
        vector<cl::Device> devices;
        devices.push_back(device);
        os << endl;
        return cl::Program(context, devices, binaries);
    } catch (cl::Error& error) {
        cerr << "Loading binary failed: " << error.what() << endl;
        exit(EXIT_FAILURE);
    }
}

cl::Kernel get_kernel(
    cl::Program& program,
    string& name)
{
    os << ">>> Loading kernel: " << name << endl;
    try {
        os << endl;
        return cl::Kernel(program, name.c_str());
    } catch (cl::Error& error) {
        cerr << "Loading kernel failed: " << error.what() << endl;
        exit(EXIT_FAILURE);
    }
}

double compute_runtime(
    cl::Event& start,
    cl::Event& end)
{
    double runtime = 0;
    runtime -= start.getProfilingInfo<CL_PROFILING_COMMAND_START>();
    runtime +=   end.getProfilingInfo<CL_PROFILING_COMMAND_START>();
    return runtime * 1e-9;
}
