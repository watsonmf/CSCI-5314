/*
* 	CSCI-5314 : HW2 : hw2.c
*	Author: Mike Watson
*
*	Reads in a FASTA file and searches through to find the most commong k-mers of a user-specified length.
*	
*	This program uses a hash table to store the number of unique k-mers found in all the sequences.
*	
*	Required libraries: POPT
*
*	Sample command line: hw1 -f example.fasta -l 5
*
*
*	Speed and Memory limitations: Overall speed should be O(nk), where n is the total length of all
*	sequences and k is the k-mer length. Since the entire file is read into memory, the program requres
*	enough memory to store all the sequences. The program will need an additional 12 * 4^k byes in order
*	to store the kmerArray hash table.
*/

#define __USE_LARGEFILE64
#define _LARGEFILE_SOURCE
#define _LARGEFILE64_SOURCE

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <popt.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <math.h>

#include "hw2.h"


int main (int argc, char** argv)
{
	char* fileName = NULL;
	int kmerLength = 0;
	
	Sequence* firstSequence = malloc(sizeof(Sequence));
	Sequence* currentSequence;
	
	/*	-------Begin POPT parsing-------	*/

	poptContext	POPT_Context;  /* context for parsing command-line options */
	int			POPT_Ret;      /* used for iterating over the arguments */

	struct poptOption optionsTable[] =
	{
		{ "file", 'f', POPT_ARG_STRING, &fileName, 'f', "Specify FASTA format genome file to read in", "FILENAME" },
		{ "length", 'l', POPT_ARG_INT, &kmerLength, 'l', "Length of K-mer to match", "LENGTH" },
		
		POPT_AUTOHELP
		{ NULL, '\0', 0, NULL, 0}
	};

	POPT_Context = poptGetContext(NULL, argc,  (const char **)argv, optionsTable, 0);
	poptSetOtherOptionHelp(POPT_Context, "[Options]\n\n[Try --help for a more detailed description of the options]\n");

	/* values are filled into the data structures by this function */
	while ((POPT_Ret = poptGetNextOpt(POPT_Context)) >= 0)
	{
		switch (POPT_Ret)
		{
		case 'f':
			/* handle file argument */
			break;
		case 'l':
			/* handle length argument */
			if (kmerLength < 3 || kmerLength > 8)
			{
				printf("Error: specified length my be between 3 and 8.\n");
				exit(0);
			}
		}
	}

	if (POPT_Ret < -1)
	{
		/* an error occurred during option processing */
		fprintf(stderr, "%s: %s\n",
		        poptBadOption(POPT_Context, POPT_BADOPTION_NOALIAS),
		        poptStrerror(POPT_Ret));
		return 1;
	}
	/*	-------End POPT parsing-------	*/
	
	if (fileName == NULL)
	{
		printf("Please specify a FASTA file for input. (%s --help to list all options)\n", argv[0]);
		
		exit(1);
	}
	
	if (kmerLength == 0)
	{
		printf("Please specify the K-mer length to search for (%s --help to list all options)\n", argv[0]);
	
		exit(1);
	}
	
	int kmerArrayLength = (int) pow(4, kmerLength);
	
	char* fileBuffer = read_file(fileName, firstSequence);
	
	Kmer* kmerArray = calloc(kmerArrayLength, sizeof(Kmer));
	
	get_kmer_count(kmerArray, firstSequence, kmerLength );

	// uses C stdlib quick sort for sort the kmerArray. 
	// Warning: index_hash function will no longer work after this sort!!
	qsort(kmerArray, kmerArrayLength, sizeof(Kmer), compare_kmer_count);
	
	int numberOfSequences = 0;
	currentSequence = firstSequence;
	
	while (currentSequence != NULL)
	{
		numberOfSequences++;
		currentSequence = currentSequence->next;
	}

	print_top_five(kmerArray, kmerArrayLength, kmerLength, numberOfSequences, fileName);
	
	return 0;
}


/*
*	Reads a file into memory and builds a linked-list containing all the sequences found in the file.
*	Uses stat64() and fopen64() to allow for loading files larger than 2GB.
*	
*	All file data is stored in the fileBuffer array and each of the sequence structs use pointers to 
*	access their respective data.
*/
char* read_file(char* inputFile, Sequence* currentSequence)
{
	char* fileBuffer;
	int fileSize;
	int c;
	int read;
	int sequenceStart;
	int index = 0;
	
	struct stat64 statStruct;
	
	// check the file size in order to malloc the necessary amount of memory to the fileBuffer array.
	if (stat64(inputFile, &statStruct) == -1) 
	{
		return NULL;
	}
	
	fileSize = statStruct.st_size;
	
	FILE* file = fopen64(inputFile, "r");

	fileBuffer = malloc(fileSize);
	
	while ((c = fgetc(file)) != '\n')
	{
		if (c == '>')
		{
			currentSequence->sequenceName = &fileBuffer[index];
		} else if (c == ' ')
		{
			fileBuffer[index++] = '\0';
			currentSequence->sequenceDescription = &fileBuffer[index];
		} else
		{
			fileBuffer[index++] = c;
		}
	}
	
	fileBuffer[index++] = '\0';
	currentSequence->sequenceType = 1;
	currentSequence->data = &fileBuffer[index];
	sequenceStart = index;
	
	// Uses fgetc to read the file in one byte at a time. 
	// Any ASCII character of lesser numerical value than 'A' will be checked if it is '>' to signify
	// a new sequence start.
	// Any other non-alphabet characters will be ignored in order to skip newlines in the sequence data.
	while((c = fgetc(file)) != EOF)
	{
		if (c < 'A')
		{
			if (c == '>')
			{
				currentSequence->sequenceLength = index - sequenceStart;
				currentSequence->next = malloc(sizeof(Sequence));
				currentSequence = currentSequence->next;
				currentSequence->sequenceName = &fileBuffer[index];
				currentSequence->sequenceType = 1;
				while ( (c = fgetc(file)) != '\n')
				{
					if (c == ' ')
					{
						fileBuffer[index++] = '\0';
						currentSequence->sequenceDescription = &fileBuffer[index];
					} else
					{
						fileBuffer[index++] = (char) c;
					}
				}
				
				fileBuffer[index++] = '\0';
				currentSequence->data = &fileBuffer[index];
				sequenceStart = index;
			} else if(c == '-')
			{
				fileBuffer[index++] = '-';
			}
		} else
		{
			fileBuffer[index++] = c;
		}
	}
	
	currentSequence->sequenceLength = index - sequenceStart;

	return fileBuffer;
}

/*
*	Checks through the linked-list of sequences and calls add_kmer for each k-mer found in the sequence.
*	Unique K-mers are added while duplicate K-mers are ignored.
*/
void get_kmer_count(Kmer* kmerArray, Sequence* currentSequence, int kmerLength )
{
	int sequenceNumber = 1;
	
	while (currentSequence != NULL)
	{
		for (int i = 0; i <= currentSequence->sequenceLength - kmerLength; i++)
		{
			add_kmer(kmerArray, &currentSequence->data[i], kmerLength, sequenceNumber);
		}
		
		currentSequence = currentSequence->next;
		sequenceNumber++;
	}
}

/*
*	Index hashing function that will read in a K-mer and return a numerical index for the kmerArray.
*	Example: AAA = 0, AAC = 1, AAG = 2, AAT = 3, ACA = 4, ... , TTT = 63 
*/
int index_hash(char* k_mer, int kmerLength)
{
	static char indexLookup[INDEX_LOOKUP_TABLE_LENGTH] = INDEX_LOOKUP_TABLE;

	int index = 0;

	for (int i = 0; i < kmerLength; i++)
	{
		if (indexLookup[k_mer[i]] < 0)
		{
			return -1;
		}
		
		index += ((int) pow(4, i)) * indexLookup[k_mer[kmerLength - 1 - i]];
	}

	return index;
}

/*
* 	Increments the count of each unique K_mer found in a sequence. Any K-mer with characters other than
*	A, C, G, T is ignored.
*/
void add_kmer(Kmer* kmerArray, char* k_mer, int kmerLength, int sequenceNumber)
{
	int index = index_hash(k_mer, kmerLength);
	
	/* ignore kmer with '-' or N or other character */
	if (index < 0)
	{
		return;
	}
	
	if (kmerArray[index].mostRecentSequence == sequenceNumber)
	{
		return;
	}
	
	kmerArray[index].mostRecentSequence = sequenceNumber;
	
	if (kmerArray[index].numberFound == 0)
	{
		kmerArray[index].k_mer = k_mer;
	}
	
	kmerArray[index].numberFound++;
}

/*
*	Prints the five most common K-mers found in all the sequences. If is a tie that extends beyond 5,
*	it will print all equally common K-mers.
*/
void print_top_five(Kmer* kmerArray, int kmerArrayLength, int kmerLength, int numberOfSequences, char* fileName)
{
	printf("5 most common %d-mer from the %d sequences in fasta file [%s]:\n", kmerLength, numberOfSequences, fileName);
	
	int i = kmerArrayLength - 1;
	
	if (kmerArray[i].k_mer == NULL)
	{
		printf("no k-mers found.\n");
		return;
	}
	
	printf("%.*s %d\n", kmerLength, kmerArray[i].k_mer, kmerArray[i].numberFound);
	
	while (--i >= kmerArrayLength - 5 || kmerArray[i].numberFound == kmerArray[i + 1].numberFound)
	{
		if (kmerArray[i].k_mer == NULL)
		{
			printf("no more k-mers found.\n");
			return;
		}
		
		printf("%.*s %d\n", kmerLength, kmerArray[i].k_mer, kmerArray[i].numberFound);
	}
}

/*
*	Compare function used for c stdlib quicksort function.
*
*	This function is used for comparison of the Kmer structs contained in the kmerArray.
*/
int compare_kmer_count(const void* kmerStruct1, const void* kmerStruct2)
{
	return (((Kmer*)kmerStruct1)->numberFound - ((Kmer*)kmerStruct2)->numberFound);
}