###########################################################################
# Compiler
###########################################################################

CC = gcc
MPICC = mpicc
CFLAGS = -std=c99 -g3
INCFLAGS = -I$(HOME)/usr/include -lm
LDFLAGS = -L$(HOME)/usr/lib -lm

ifeq ($(DEBUG),on)
	CFLAGS += -DDEBUG
endif

###########################################################################
# File Declarations
###########################################################################
BINARIES = hw4

# C, H, and Object Files
C_FILES = watson_hw4.c
H_FILES = watson_hw4.h
O_FILES = watson_hw4.o

MPI_C_FILES =
MPI_H_FILES =
MPI_O_FILES =

###########################################################################
# Main Targets
###########################################################################
all: $(BINARIES)

janus: $(BINARIES)

hw4: $(O_FILES)
	$(CC) $(INCFLAGS) -o ./hw4 $(O_FILES) $(LDFLAGS)

clean:
	/bin/rm -rf core $(O_FILES) $(MPI_O_FILES) $(BINARIES)

# Suffixes

.c.o:
	$(CC) $(CFLAGS) $(INCFLAGS) -c -o $*.o $*.c $(LDFLAGS)

# Dependencies
