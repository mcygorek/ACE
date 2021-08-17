
ifneq ($(MKLROOT),)
  MKL_FLAGS=-DMKL_ILP64 -m64 -I$(MKLROOT)/include -L$(MKLROOT)/lib/intel64 -Wl,--no-as-needed -lmkl_intel_ilp64 -lmkl_gnu_thread -lmkl_core -lgomp -lpthread -lm -ldl -DEIGEN_USE_MKL_ALL
#  $(info $$MKL_FLAGS is [${MKL_FLAGS}])
else 
  $(info )
  $(warning ACE will not be compiled against MKL. If you would like to use it, please specify variable MKLROOT, e.g., MKLROOT=/opt/intel/mkl/. Using MKL routines for SVDs is often faster but not as accurate, which may lead to unphysical results for very small thresholds.)
  $(info )
endif

ifeq ($(EIGEN_HOME),)
  EIGEN_HOME=/usr/include/eigen3/
endif
ifeq (,$(wildcard $(EIGEN_HOME)/Eigen/Eigen))
    $(error Please specify variable EIGEN_HOME)
endif

#FFTW_FLAGS=-lfftw3

OPTS += -O3 $(FFTW_FLAGS) -lm -g -m64


