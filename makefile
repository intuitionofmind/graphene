INC = -I/vol-th/home/lattice/zhengwei/fftw/include/ -I/usr/local/mpi-gcc46/include/
LIB = -lm -L/vol-th/home/lattice/zhengwei/fftw/lib/ -lfftw3 -L/usr/local/mpi-gcc46/lib/

CXX = mpicxx 
CXXFLAGS = $(INC) -O2 -Wall -ansi -funroll-loops -std=gnu++0x 
OBJS = main.o operation.o extra.o distribution.o twopoint.o ##acceleration.o check.o spin_cor.o

TARGET = graphene

all:$(TARGET)
$(TARGET):$(OBJS)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LIB)

main.o:main.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@ $(LIB)
operation.o:operation.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@ $(LIB)
extra.o:extra.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@ $(LIB) 
distribution.o:distribution.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@ $(LIB) 
twopoint.o:twopoint.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@ $(LIB)
#spin_cor.o:spin_cor.cpp
#	$(CXX) $(CXXFLAGS) -c $< -o $@ $(LIB)
#acceleration.o:acceleration.cpp
#	$(CXX) $(CXXFLAGS) -c $< -o $@ $(LIB)
#check.o:check.cpp
#	$(CXX) $(CXXFLAGS) -c $< -o $@ $(LIB)

clean:
	rm -f $(TARGET) *.o *.log sample *.dat
