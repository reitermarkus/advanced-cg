TARGET = Radiosity
TARGET_TEST = test/run

OBJ = $(filter-out $(TARGET).o, $(patsubst %.cpp, %.o, $(wildcard *.cpp)))
OBJ_TEST = $(patsubst %.cpp, %.o, $(wildcard test/*.cpp))

INCLUDES = -Iinclude -I./

CXXFLAGS += -O3 -Wall -std=c++1z -fopenmp $(INCLUDES)

ifeq ($(shell command -v module 2>/dev/null), module)
  CXXFLAGS += -static-libstdc++
  CXX := module load gcc/7.2.0 || true && $(CXX)
endif

# Rules
all: $(TARGET)

$(TARGET): $(TARGET).o $(OBJ)
	$(CXX) $(CXXFLAGS) $^ -o $@

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $^ -o $@

.PHONY: test
test: $(OBJ_TEST) $(OBJ)
	$(CXX) $(CXXFLAGS) $^ -o $(TARGET_TEST)
	./$(TARGET_TEST)

.PHONY: clean
clean:
	$(RM) $(TARGET) *.o */*.o

run: all
	./$(TARGET)
