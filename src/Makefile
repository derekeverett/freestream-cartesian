#compiler
CC = g++

#compiler flags
# -g adds debug info
# -Wall turns on most warnings
CFLAGS = -g -Wall
LIBS= -lgslcblas -lgsl
#build target
TARGET = Run

all: $(TARGET)

$(TARGET): $(TARGET).cpp
	$(CC) $(CFLAGS) $(LIBS) -o $(TARGET) $(TARGET).cpp
clean:
