OBJ = $(patsubst %.cpp, %.o, $(wildcard *.cpp))
TARGET = Radiosity

CXXFLAGS = -O3 -Wall -std=c++1z

ifeq ($(shell command -v module 2> /dev/null), module)
  CXX := module load gcc/7.2.0 || true && ${CXX}
endif

# Rules
all: $(TARGET)

$(TARGET): $(OBJ)
	$(CXX) $(CXXFLAGS) $^ -o $@

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $^ -o $@

.PHONY: clean
clean:
	$(RM) $(TARGET) *.o

run: all
	./$(TARGET)
