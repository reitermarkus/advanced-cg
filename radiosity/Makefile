OBJ = $(patsubst %.cpp, %.o, $(wildcard *.cpp))
TARGET = Radiosity

CXXFLAGS = -O3 -Wall -std=c++1y

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
