/*
* 	CSCI-5314 : HW1 : hw1.c
*	Author: Mike Watson
*
*	Uses the following binary bytes for reducing conditionals within the search loop using binary AND:
*
*	A 00000001
*	C 00000010
*	G 00000100
*	T 00001000
*	B 00001110
*	D 00001101
*	H 00001011
*	K 00001100
*	M 00000011
*	N 00001111
*	R 00000101
*	S 00000110
*	V 00000111
*	W 00001001
*	Y 00001010
*	
*	Protein sequences are matched directly withough XAND operations.
*
*/

#define LOOKUP_TABLE_LENGTH 124
#define LOOKUP_TABLE {0,0,0,0,0,0,0,0,0,0, 0 /*/n*/,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,   \
0 ,0, 0, 0, 0, 0, 0, 0, 0, 0, 0 , 0, 0 /* - */, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /* > */, 0, 0, \
1, 46, 2, 13, 65, 66, 4, 11, 67, 97, 12, 68, 3, 32, 69, 70, 71, 5, 6, 8, 24, 7, 9, 98, 10, 99, 0, 0, 0, 0, 0, 0,  \
1, 46, 2, 13, 65, 66, 4, 11, 67, 97, 12, 68, 3, 32, 69, 70, 71, 5, 6, 8, 24, 7, 9, 98, 10, 99} 

#define KMER_LOOKUP_TABLE_LENGTH 124
#define KMER_LOOKUP_TABLE {0,0,0,0,0,0,0,0,0,0, 0 /*/n*/,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,   \
0 ,0 ,0, 0, 0, 0, 0, 0, 0, 0, 0 , 0, 0 /* - */, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /* > */, 0, 0, \
1, 14, 2, 13, 0, 0, 4, 11, 0, 0, 12, 0, 3, 15, 0, 0, 0, 5, 6, 8, 8, 7, 9, 0, 10, 0, 0, 0, 0, 0, 0, 0,  \
1, 14, 2, 13, 0, 0, 4, 11, 0, 0, 12, 0, 3, 15, 0, 0, 0, 5, 6, 8, 8, 7, 9, 0, 10, 0} 

#define REVERSE_COMPLIMENT_TABLE_LENGTH 16
#define REVERSE_COMPLIMENT_TABLE {0, 8, 4, 12, 2, 10, 6, 14, 1, 9, 5, 13, 3, 11, 17, 15}

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

char* read_file(char* inputFile, Sequence* currentSequence, char* lookupTable);
void search_nucleic_acid(Sequence* currentSequence, char* k_mer, char* reverseComplementKmer, int k);
void search_amino_acid(Sequence* currentSequence, char* k_mer, int k);
void print_match(Sequence* currentSequence, int start, int end, char strand);

