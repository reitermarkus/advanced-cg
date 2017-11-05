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

$(TARGET): $(OBJ)

$(TARGET).o: $(TARGET).cpp
	$(CC) $(CFLAGS) $(INCLUDES) -c $^ -o $@

%.o: %.cpp
	$(CC) $(CFLAGS) $(INCLUDES) -c $^ -o $@

.PHONY: clean
clean:
	$(RM) $(TARGET) *.o

run: all
	./$(TARGET)


