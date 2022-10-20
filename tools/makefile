NETCDF = -lnetcdf_c++ -lnetcdf
SRC_FILES := $(wildcard ../src/*.cpp)
CXX = g++
CXXFLAGS = -Wall -ggdb

forefire2: $(SRC_FILES)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(NETCDF)