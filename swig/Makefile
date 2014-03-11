# SWIG
INTERFACES = PLibForeFire.i
WRAPPERS   = $(INTERFACES:.i=_wrap.cxx)
PROXIES    = $(INTERFACES:.i=.py      )

# Default target: build the tests
.PHONY : all
all: $(WRAPPERS) PLibForeFire.cpp PLibForeFire.h
#	CC=g++ ./setup.py build_ext -i
	g++ -I/Users/filippi/soft/lib/python/numpy/core/include -I/Users/filippi/soft/include/ -I/System/Library/Frameworks/Python.framework/Versions/2.7/include/python2.7 -fmessage-length=0  -c PLibForeFire.cpp PLibForeFire_wrap.cxx  
	g++ -L/Users/filippi/soft/lib -lpython  -shared -lnetcdf -lnetcdf_c++ -lForeFire  PLibForeFire.o PLibForeFire_wrap.o -o _forefire.so
.PHONY : tense
tense: $(WRAPPERS) Tensor.cxx Tensor.h
	./setuptense.py build_ext -i


# Rule: %.i -> %_wrap.cxx
%_wrap.cxx: %.i %.h numpy.i
	swig -c++ -python $<
%_wrap.cxx: %.i %1.h %2.h numpy.i
	swig -c++ -python $<

# Clean target
.PHONY : clean
clean:
	$(RM) -r build
	$(RM) *.so
	$(RM) *.o
	$(RM) forefire.py*
	$(RM) $(WRAPPERS)
	$(RM) $(PROXIES)
