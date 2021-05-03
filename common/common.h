#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>

#define CL_HPP_ENABLE_EXCEPTIONS
#define CL_HPP_MINIMUM_OPENCL_VERSION 120
#define CL_HPP_TARGET_OPENCL_VERSION 120
#define CL_HPP_ENABLE_PROGRAM_CONSTRUCTION_FROM_ARRAY_COMPATIBILITY
#include <CL/cl2.hpp>

void init(
    cl::Context &context,
    std::vector<cl::Device> &devices);

void print_platform(
    cl::Platform &platform);

void print_device(
    cl::Device &device,
    bool marker = false);

std::string get_source(
    std::string& filename);

std::string get_flags();

cl::Program compile_program(
    cl::Context& context,
    cl::Device& device,
    std::string& source);

void write_source(
    std::string& source,
    std::string& filename);

cl::Program get_program(
    cl::Context& context,
    cl::Device& device,
    std::string& filename);

cl::Kernel get_kernel(
    cl::Program& program,
    std::string& name);

double compute_runtime(
    cl::Event& start,
    cl::Event& end);
