
#define INDEX_LOOKUP_TABLE_LENGTH 124
#define INDEX_LOOKUP_TABLE {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, -1 /*/n*/,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,   \
-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1 /* - */,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1 /* > */,-1,-1, \
0,-1,1,-1,-1,-1,2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,3,3,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,  \
0,-1,1,-1,-1,-1,2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,3,3,-1,-1,-1,-1,-1}

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

typedef struct _KMER
{
	int numberFound;
	char* k_mer;
	int mostRecentSequence;
} Kmer;


char* read_file(char* inputFile, Sequence* currentSequence);
void get_kmer_count(Kmer* kmerArray, Sequence* currentSequence, int kmerLength );
int index_hash(char* k_mer, int kmerLength);
void add_kmer(Kmer* kmerArray, char* k_mer, int kmerLength, int sequenceNumber);
void print_top_five(Kmer* kmerArray, int kmerArrayLength, int kmerLength, int numberOfSequences, char* fileName);
void test_print(Kmer* kmerArray, int kmerArrayLength, int kmerLength);
void test_print_all(char* fileBuffer);
int compare_kmer_count(const void* kmerStruct1, const void* kmerStruct2);
