#compiler
CC = g++-7

#compiler flags
# -g adds debug info
# -Wall turns on most warnings
CFLAGS = -g -Wall -fopenmp #comment or remove openmp if parallelization is not possible
LIBS= -lgslcblas -lgsl
#build target
SRC = src
TARGET = Run

all: $(TARGET)

$(TARGET): $(SRC)/$(TARGET).cpp
	$(CC) $(CFLAGS) $(LIBS) -o $(TARGET) $(SRC)/$(TARGET).cpp
clean:
