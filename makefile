INC = -I/home/skysniper/fftw/include/
LIB = -lm -L/home/skysniper/fftw/lib/ -lfftw3

CXX = g++
CXXFLAGS = $(INC) -fopenmp -O2 -Wall -pedantic -ansi #-std=c++11 
OBJS = main.o function_class.o distribution.o twopoint.o ##acceleration.o check.o spin_cor.o

TARGET = Npeak2

all:$(TARGET)
$(TARGET):$(OBJS)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LIB)

main.o:main.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@ $(LIB)
function_class.o:function_class.cpp
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
