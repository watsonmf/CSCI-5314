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
	int matrixWidth;
	int matrixHeight;
} AlignmentMatrix;

typedef struct _ScoringMatrix
{
	int numChars;
	int* matrix;
	char* lookupTable;
	char* lookupTableReverse;
} ScoringMatrix;

typedef struct _Node
{
	struct _Node* left;
	struct _Node* right;
	float distance;
	float height;
	Sequence* sequence;
} Node;


ScoringMatrix* read_scoring_file(char* scoreFile);
Node* build_tree(ScoringMatrix* scoreMatrix, Sequence* currentSequence, int gapPenalty, int numberOfSequences);
char* read_file(char* inputFile, Sequence* currentSequence, ScoringMatrix* score, int* numberOfSequences);
float global_alignment_distance(Sequence* sequenceOne, Sequence* sequenceTwo, ScoringMatrix* score, int gapPenalty);
void print_matrix(AlignmentMatrix* alignment);
void print_tree(Node* aNode);
void print_tree_recursive(Node* aNode);

