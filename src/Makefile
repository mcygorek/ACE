
ifneq ($(MKLROOT),)
  MKL_FLAGS=-DMKL_ILP64 -m64 -I$(MKLROOT)/include -L$(MKLROOT)/lib/intel64 -Wl,--no-as-needed -lmkl_intel_ilp64 -lmkl_gnu_thread -lmkl_core -lgomp -lpthread -lm -ldl -DEIGEN_USE_MKL_ALL
#  $(info $$MKL_FLAGS is [${MKL_FLAGS}])
else 
  $(info )
  $(warning Please specify variable MKLROOT, e.g., MKLROOT=/opt/intel/mkl/. Not using MKL will lead to very slow calculation of SVDs!)
  $(info )
endif

ifeq ($(EIGEN_HOME),)
  EIGEN_HOME=/usr/include/eigen3/
endif
ifeq (,$(wildcard $(EIGEN_HOME)/Eigen/Eigen))
    $(error Please specify variable EIGEN_HOME)
endif

#FFTW_FLAGS=-lfftw3

OPTS += -O3 $(FFTW_FLAGS) -lm -DPRINT_MAX_DIM

OBJS = ACE ACE_integrate_t1 
OBJS += iQUAPI TEMPO 
OBJS += readexpression
OBJS += PT_extend spectra_from_dynamics

all: $(OBJS)


%: %.cpp
	$(CXX) -o $@ $< -I../include -I$(EIGEN_HOME) $(OPTS) $(MKL_FLAGS)
	mkdir -p ../bin && mv $@ ../bin/


.PHONY: clean

clean:
	rm results.dat $(OBJS) 2>/dev/null || true
