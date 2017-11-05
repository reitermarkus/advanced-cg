CC = c++
LD = c++

OBJ = $(patsubst %.cpp, %.o, $(wildcard *.cpp))
TARGET = Radiosity

CFLAGS = -O3 -Wall -std=c++1y
LDLIBS =
INCLUDES =

SRC_DIR =
BUILD_DIR =
VPATH =

# Rules
all: $(TARGET)

$(TARGET).o: $(TARGET).cpp
	$(CC) $(CFLAGS) $(INCLUDES) -c $^ -o $@

%.o: %.cpp
	$(CC) $(CFLAGS) $(INCLUDES) -c $^ -o $@

clean:
	rm -f *.o $(TARGET)

.PHONY: clean

# Dependencies
$(TARGET): $(OBJ)

