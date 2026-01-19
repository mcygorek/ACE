$(info Compiler is $(CXX))
$(info Operating system is $(shell uname -o))

ifneq ($(MKLROOT),)
  MKL_INCL = -I$(MKLROOT)/include -DEIGEN_USE_MKL_ALL
  MKL_LINK = -DMKL_LP64 -m64 -L$(MKLROOT)/lib/intel64 -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core -lgomp -lpthread -lm -ldl 
else 
  $(info )
  $(info ########################################)
  $(info ACE will not be compiled against MKL. If you would like to use it, please specify variable MKLROOT, e.g., MKLROOT=/opt/intel/mkl/. )
  $(info ########################################)
  $(info )
endif


ifeq ($(EIGEN_HOME),)
  EIGEN_HOME=/usr/include/eigen3/
endif
ifeq (,$(wildcard $(EIGEN_HOME)/Eigen/Eigen))
    $(error Please specify variable EIGEN_HOME)
endif

#FFTW_FLAGS=-lfftw3

ifeq ($(OPTIM),)
  OPTIM = 3
endif

LIBNAME=libACE.so
ifeq ($(shell uname -o),Msys) 	
  OPTS += -DM_PI=3.14159265358979323846
  LIBNAME=libACE.dll
endif

OPTS += -O$(OPTIM) -fno-math-errno -g -m64 --std=c++11
ifneq ($(CXX),icpc)
  OPTS += -Wno-ignored-attributes
endif
OPTS += -Winvalid-pch -fPIC -pthread
OPTS += -Iinclude -I$(EIGEN_HOME) $(MKL_INCL)
OPTSLINK += $(OPTS) -lm $(MKL_LINK)

LIBSRC = $(wildcard src/*.cpp)
LIBOBJS = $(patsubst src/%.cpp, lib/%.o, $(LIBSRC))
LIBSTRING = '$(shell pwd)/lib'

EXECOBJS = ACE QUAPI TEMPO #ACE_Network #ACE_env_obs 
BINEXEC = $(patsubst %, bin/%, $(EXECOBJS))

PYBINDSUF = $(shell python3-config --extension-suffix)
PYBINDINC = $(shell python3 -m pybind11 --includes)
#PYBINDSUF = $(shell python3 -c "import sysconfig; print(sysconfig.get_config_var('EXT_SUFFIX'))" )

#$(info $(PYBINDSUF) $(PYBINDINC) )

all: code tools
code: PCH $(EXECOBJS) 
$(EXECOBJS): %: src_exec/%.cpp  lib/$(LIBNAME) 
	$(CXX) -o bin/$@ $< $(OPTSLINK) -Llib -lACE -Wl,-rpath,$(LIBSTRING)

pybind: pybind/ACE$(PYBINDSUF)

pybind/ACE$(PYBINDSUF): pybind/pybind.cpp lib/$(LIBNAME)
	$(CXX) -O3 -shared -std=c++11 -fPIC  -Iinclude -I$(EIGEN_HOME) $(MKL_INCL) $(PYBINDINC) pybind/pybind.cpp -o pybind/ACE$(PYBINDSUF) $(OPTSLINK) -Llib -lACE -Wl,-rpath,$(LIBSTRING)
#	$(CXX) -O3 -Wall -shared -std=c++11 -fPIC $(PYBINDINC) pybind/pybind.cpp -o pybind/ACE$(PYBINDSUF) $(OPTSLINK) -Llib -lACE -Wl,-rpath,$(LIBSTRING) -Iinclude -I$(EIGEN_HOME) $(MKL_INCL)

lib: lib/$(LIBNAME)
lib/$(LIBNAME): $(LIBOBJS) 
	$(CXX) -shared -fPIC -o lib/$(LIBNAME) -flto lib/*.o 

$(LIBOBJS): lib/%.o: src/%.cpp include/%.hpp 
	$(CXX) -o $@ -c $< -fPIC $(OPTS)

PCH: include/PCH.hpp.gch 
include/PCH.hpp.gch: include/PCH.hpp
	$(CXX) -x c++-header -o include/PCH.hpp.gch -c include/PCH.hpp $(OPTS)

#TOOLS = PT_build PT_sweep PT_join 
#TOOLS += PT_analyze PT_apply_system_propagator PT_coarse_grain PT_extend
TOOLS = integrate_outfile differentiate_outfile FT_outfile smoothen_outfile
TOOLS += PTB_sweep_forward PTB_sweep_backward PTB_apply_system_propagator
TOOLS += PTB_analyze PTB_join PTR_to_PTB PTB_strip_env_ops PTB_extend
TOOLS += PTB_transfer_env_ops PTB_expand_in_place PTB_coarse_grain 
TOOLS += PTB_outer_contract
TOOLS += PTR_sweep_backward PTR_analyze
TOOLS += PME_renorm BCF_from_J
TOOLS += readexpression print_densmat timedep_eigenstates

EXPERIMENTAL = generate_J generate_noise 
EXPERIMENTAL += DiagBB_print_K single_mode_K block_combine 
EXPERIMENTAL += extract_singular_values 
EXPERIMENTAL += richardson_extrapolate #Chebyshev_Expand
EXPERIMENTAL += TC_overlap TC_PT TC_join TC_TinvTcombine test_CompressionTree
EXPERIMENTAL += test_FFT estimate_memory test_split
EXPERIMENTAL += test_buffer test_GaussNewton fit_K_single_mode 
EXPERIMENTAL += test_MeierTannor test_DrudeLorentz
EXPERIMENTAL += extract_effective_propagator DynamicalMap PT_traceout
EXPERIMENTAL += test_HermitianLiouvilleBasis test_Largest_EV test_Largest_EV2
EXPERIMENTAL += test_QUAPI_vs_PT test_RandomizedCompression

tools: $(TOOLS)
$(TOOLS): %: src_exec/tools/%.cpp lib/$(LIBNAME)
	$(CXX) -o tools/$@ $< $(OPTSLINK) -Llib -lACE -Wl,-rpath,$(LIBSTRING)

experimental: $(EXPERIMENTAL)
$(EXPERIMENTAL): %: src_exec/experimental/%.cpp lib/$(LIBNAME)
	$(CXX) -o tools/$@ $< $(OPTSLINK) -Llib -lACE -Wl,-rpath,$(LIBSTRING)


.PHONY: clean

clean:
	rm -rf bin/* tools/* lib/* include/PCH.hpp.gch pybind/*.so

