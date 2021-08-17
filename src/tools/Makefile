include ../../make_params.mk

OBJS = generate_J
OBJS += DiagBB_print_K  
OBJS += PT_analyze PT_apply_system_propagator PT_coarse_grain PT_extend

all: $(OBJS)


%: %.cpp
	$(CXX) -o $@ $< -I../../include -I$(EIGEN_HOME) $(OPTS) $(MKL_FLAGS)
	mkdir -p ../../tools/ && mv $@ ../../tools/


.PHONY: clean

clean:
	rm results.dat $(OBJS) 2>/dev/null || true

