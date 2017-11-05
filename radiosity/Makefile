CC = c++
LD = c++

OBJ = Radiosity.o Vector.o Image.o Rectangle.o
TARGET = Radiosity

CFLAGS = -O3 -Wall -std=c++11
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

