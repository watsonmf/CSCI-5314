/*
* 	CSCI-5314 : HW3 : watson_hw3.h
*	Author: Mike Watson
*

/* definitions for global alignment penalties */
#define GLOBAL_GAP 1
#define GLOBAL_MISMATCH 1
#define GLOBAL_MATCH 0

/* definitions for local alignment penalties */
#define LOCAL_GAP -1
#define LOCAL_MISMATCH -1
#define LOCAL_MATCH 3

/* definitions for direction values used in 2d matrix */
#define LEFT	0
#define DIAG	1
#define UP		2

typedef char type;

typedef struct _SEQUENCE
{
	struct _SEQUENCE* next;
	char* sequenceName;
	char* sequenceDescription;
	int sequenceLength;
	type sequenceType; // 1 = nucleic acid, 2 = amino acid
	char* data; // pointer to the start of sequence.
} Sequence;

typedef struct _AlignmentMatrix
{
	int* score;
	char* direction;
	int matrixWidth;
	int matrixHeight;
} AlignmentMatrix;


char* read_file(char* inputFile, Sequence* currentSequence);
AlignmentMatrix* global_alignment(Sequence* sequenceOne, Sequence* sequenceTwo);
AlignmentMatrix* local_alignment( Sequence* sequenceOne, Sequence* sequenceTwo );
void print_global_alignment(AlignmentMatrix* alignment, Sequence* sequenceOne, Sequence* sequenceTwo);
void print_local_alignment(AlignmentMatrix* alignment, Sequence* sequenceOne, Sequence* sequenceTwo);
void print_matrix(AlignmentMatrix* alignment);

