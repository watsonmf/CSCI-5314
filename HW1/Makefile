########################################################################### 
# Compiler
########################################################################### 

CC = gcc
MPICC = mpicc
CFLAGS = -std=c99
INCFLAGS = -I$(HOME)/usr/include -lpopt
LDFLAGS = -L$(HOME)/usr/lib -lpopt -fopenmp

ifeq ($(DEBUG),on)
	CFLAGS += -DDEBUG
endif

########################################################################### 
# File Declarations
########################################################################### 
BINARIES = hw1

# C, H, and Object Files
C_FILES = hw1.c
H_FILES = 
O_FILES = hw1.o

MPI_C_FILES = 
MPI_H_FILES = 
MPI_O_FILES = 

########################################################################### 
# Main Targets
########################################################################### 
all: $(BINARIES)

janus: $(BINARIES)

hw1-serial: $(O_FILES)
	$(CC) $(INCFLAGS) -o ./bin/hw1 $(O_FILES) $(LDFLAGS)

clean:
	/bin/rm -rf core $(O_FILES) $(MPI_O_FILES) $(BINARIES)

# Suffixes

.c.o:		
	$(MPICC) $(CFLAGS) $(INCFLAGS) -c -o $*.o $*.c $(LDFLAGS)

# Dependencies
