# Compiler flags
CXX = g++
CXXFLAGS = -g -fopenmp -std=c++11 -O2 -I..
LDFLAGS = -g -fopenmp -std=c++11

#INCLUDES += -I${OPENCL_INCLUDE}
#LDFLAGS += -lOpenCL -L${OPENCL_LIB}
#INCLUDES += -I/cm/shared/package/cuda101/toolkit/10.1.105/targets/x86_64-linux/include/
#LDFLAGS += -L/cm/shared/package/cuda101/toolkit/10.1.105/targets/x86_64-linux/lib -lOpenCL

# Altera OpenCL compiler flags
INCLUDES += $(shell aocl compile-config)
LDFLAGS += $(shell aocl link-config)
#LDFLAGS += -L/usr/local/intel/fpga/inteldevstack/intelFPGA_pro/hld/linux64/lib/ -lalteracl
AOC = aoc

ifneq ("$(USE_SIN_COS_LOOKUP_TABLE)", "")
AOCOFLAGS += -DUSE_SIN_COS_LOOKUP_TABLE
endif

ifneq ("$(NR_GRIDDERS)", "")
AOCOFLAGS += -DNR_GRIDDERS=$(NR_GRIDDERS)
endif

ifneq ("$(SEED)", "")
AOCXFLAGS += -seed=$(SEED)
endif

ifneq ("$(PARALLEL)", "")
AOCRFLAGS += -parallel=$(PARALLEL)
else
AOCRFLAGS += -parallel=6
endif

#AOCOFLAGS +=-board=p385a_min_ax115 -DAOCL_BOARD_p385a_min_ax115
#AOCOFLAGS += -board=pac_s10_dc -DAOCL_BOARD_pac_s10_dc
AOCOFLAGS += -I${INTELFPGAOCLSDKROOT}/include/kernel_headers
#AOCOFLAGS += -ffp-reassoc -ffp-contract=fast
#AOCOFLAGS += -I/home/romein/projects/Triple-A/FPGA/IDG/fpga
#AOCRFLAGS += -L/home/romein/projects/Triple-A/FPGA/IDG/fpga/fft_lib -lifft_32x32_s10.aoclib
AOCOFLAGS += -cl-fast-relaxed-math
AOCRFLAGS += -global-ring
AOCRFLAGS += -no-hardware-kernel-invocation-queue
#AOCOFLAGA += -DHYPER_OPTIMIZED_HANDSHAKING_DISABLED
#AOCRFLAGS += -hyper-optimized-handshaking=off
#AOCRFLAGS += -ecc
AOCRFLAGS += -report
#AOCRFLAGS += -profile=all
AOCXFLAGS += ${SEEDFLAGS}
AOCRFLAGS += -clock=400MHz
AOCXFLAGS += -high-effort

# Emulator
#AOCOFLAGS += -march=emulator -legacy-emulator -DEMULATOR
#AOCRFLAGS += -emulator-channel-depth-model=strict
# export CL_CONTEXT_EMULATOR_DEVICE_INTELFPGA=1

# FFTW flags
FFTW_LDFLAGS = -lfftw3f

CLSOURCES = degridder.cl fft.cl gridder.cl sincos.cl
CXXSOURCES = common/common.cpp common/init.cpp reference/degridder.cpp reference/fft.cpp reference/gridder.cpp run-fpga-degridder.cpp run-fpga-gridder.cpp run-fpga-sincos.cpp run-gpu-degridder.cpp run-gpu-gridder.cpp run-fpga-fft.cpp
OBJECT_FILES = $(CXXSOURCES:%.cpp=%.o)
EXECUTABLES = $(CXXSOURCES:%.cpp=%.x)
FPGA_STUFF_SUBDIRS = $(CLSOURCES:%.cl=%)
DEPENDENCIES = $(CXXSOURCES:%.cpp=%.d) $(CLSOURCES:%.cl=%.d)

TMPDIR ?= /tmp

default:: run-fpga-gridder.x gridder.aocx run-fpga-degridder.x degridder.aocx run-fpga-fft.x fft.aocx run-fpga-sincos.x sincos.aocx

run-fpga-%.x: run-fpga-%.o common/common.o common/init.o reference/fft.o reference/%.o
	${CXX} ${LDFLAGS} ${FFTW_LDFLAGS} -o $@ $^

run-fpga-fft.x: run-fpga-fft.o common/common.o reference/fft.o
	${CXX} ${LDFLAGS} ${FFTW_LDFLAGS} -o $@ $^

run-fpga-sincos.x: run-fpga-sincos.o common/common.o
	${CXX} ${LDFLAGS} -o $@ $^

%.d: %.cpp
	${CXX} -MM -MT $*.o ${CXXFLAGS} ${INCLUDES} -o $@ $^

%.d: %.cl
	${CXX} -MM -MT $*.aoco -x c ${INCLUDES} -o $@ $^

%.o: %.cpp
	${CXX} -c ${CXXFLAGS} ${INCLUDES} -o $@ $<

%.aoco: $(SRCDIR)/%.cl
	rm -rf $*
	${AOC} -c ${AOCOFLAGS} -o $@ $<

%.aoco: %.cl
	rm -rf $*
	${AOC} -c ${AOCOFLAGS} -o $@ $<

%.aocr: %.aoco
	${AOC} -rtl ${AOCRFLAGS} -o $@ $<

%.aocx: %.aocr
	(unset DISPLAY; ${AOC} ${AOCXFLAGS} -o $@ $<)

%/build:
	test -f $@ || test -f /tmp/stop || (mkdir -p $* && (echo `hostname` && mkdir -p /local-ssd/`basename $*` && cd /local-ssd/`basename $*` && SRCDIR=$(PWD) $(if $(findstring LU, $(word 2, $(subst _, , $*))), USE_SIN_COS_LOOKUP_TABLE=1) NR_GRIDDERS=$(word 3, $(subst _, , $*)) SEED=$(word 4, $(subst _, , $*)) time make -f $(PWD)/Makefile $(word 1, $(subst _, ,$(notdir $*))).aoco $(word 1, $(subst _, ,$(notdir $*))).aocr $(word 1, $(subst _, ,$(notdir $*))).aocx && awk '/ MHz/ || /freq:/' $(word 1, $(subst _, ,$(notdir $*)))/quartus_sh_compile.log|tail -n 1; cd $(PWD); mv /local-ssd/`basename $*`/* $*; rmdir /local-ssd/`basename $*` )>$@ 2>&1)

clean::
	rm -rf $(EXECUTABLES) $(OBJECT_FILES) $(FPGA_STUFF_SUBDIRS) $(DEPENDENCIES)

sleep.%::
	$(subst ., ,$@)

ifeq (0, $(words $(findstring $(MAKECMDGOALS), clean)))
-include $(DEPENDENCIES)
endif
