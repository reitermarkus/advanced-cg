TARGET = Corona

OBJ = $(filter-out $(TARGET).o, $(patsubst %.cpp, %.o, $(wildcard *.cpp)))
OBJ_SHARED = $(patsubst ../shared/%.cpp, ../shared/%.o, $(wildcard ../shared/*.cpp))

INCLUDES = -I./ -I../shared
CXXFLAGS = -O3 -Wall -Wextra -std=c++1z -fopenmp


# Rules
all: $(TARGET)

$(TARGET): $(TARGET).o $(OBJ) $(OBJ_SHARED)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $^ -o $@

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $^ -o $@

.PHONY: clean
clean:
	$(RM) $(TARGET) $(TARGET_TEST) *.o */*.o ../shared/*.o

run: all
	time ./$(TARGET)
