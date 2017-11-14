TARGET = Radiosity
TARGET_TEST = test/run

OBJ = $(filter-out $(TARGET).o, $(patsubst %.cpp, %.o, $(wildcard *.cpp)))
OBJ_TEST = $(patsubst %.cpp, %.o, $(wildcard test/*.cpp))

INCLUDES = -Iinclude -I./
CXXFLAGS = -O3 -Wall -Wextra -std=c++1z -fopenmp

ifneq ($(OS), Windows_NT)
	ifeq ($(shell uname -s), Darwin)
		ifeq ($(shell brew ls llvm &>/dev/null; echo $$?), 0)
			LLVM_PREFIX = $(shell brew --prefix llvm)
			CXX = $(LLVM_PREFIX)/bin/clang++
			INCLUDES += -I$(LLVM_PREFIX)/include
			LDFLAGS += -L$(LLVM_PREFIX)/lib
		else
			$(error "LLVM is not installed, run `brew install llvm`")
		endif
	else
		ifeq ($(shell module avail gcc 2>&1 | grep -q gcc/7.2.0; echo $$?), 0)
			CXXFLAGS += -static-libstdc++
			CXX := module load gcc/7.2.0 && $(CXX)
		endif
	endif
else
	CXX = g++-7.2.0
endif

ifeq ($(shell ccache -V &>/dev/null; echo $$?), 0)
	CXX := ccache $(CXX)
endif

# Rules
all: $(TARGET)

$(TARGET): $(TARGET).o $(OBJ)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $^ -o $@

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $^ -o $@

.PHONY: test
test: $(OBJ_TEST) $(OBJ)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $^ -o $(TARGET_TEST)
	./$(TARGET_TEST)

.PHONY: clean
clean:
	$(RM) $(TARGET) *.o */*.o

run: all
	time ./$(TARGET)
