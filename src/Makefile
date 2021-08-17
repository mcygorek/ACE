include ../make_params.mk


OBJS = ACE iQUAPI TEMPO 

all: $(OBJS)


%: %.cpp
	$(CXX) -o $@ $< -I../include -I$(EIGEN_HOME) $(OPTS) $(MKL_FLAGS)
	mkdir -p ../bin && mv $@ ../bin/


.PHONY: clean

clean:
	rm results.dat $(OBJS) 2>/dev/null || true
