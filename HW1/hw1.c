/*
* 	CSCI-5314 : HW1 : hw1.c
*	Author: Mike Watson
*
*	Reads in a FASTA file and searches through to match the given K-mer.
*
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

#include "hw1.h"

int main (int argc, char** argv)
{
	char* fileName = NULL;
	char* chromosome = NULL;
	char* k_mer = NULL;
	Sequence* firstSequence = malloc(sizeof(Sequence));
	Sequence* currentSequence;
	
	char lookupTable[LOOKUP_TABLE_LENGTH] = LOOKUP_TABLE ;
	char k_merLookupTable[KMER_LOOKUP_TABLE_LENGTH] = KMER_LOOKUP_TABLE ;
	char reverseComplementTable[REVERSE_COMPLIMENT_TABLE_LENGTH] = REVERSE_COMPLIMENT_TABLE ;
	/*	-------Begin POPT parsing-------	*/

	poptContext	POPT_Context;  /* context for parsing command-line options */
	int			POPT_Ret;      /* used for iterating over the arguments */

	struct poptOption optionsTable[] =
	{
		{ "file", 'f', POPT_ARG_STRING, &fileName, 'f', "Specify FASTA format genome file to read in", "FILENAME" },
		{ "chromosome", 'c', POPT_ARG_STRING, &chromosome, 'c', "Specify sequence to search (default is all)", "NAME" },
		{ "k-mer", 'k', POPT_ARG_STRING, &k_mer, 'k', "K-mer to search genome for", "K-MER" },
		
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
		case 'c':
			/* handle chromosome argument */
			
			break;
		case 'k':
			/* handle k-mer argument */
			break;
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
	
	if (k_mer == NULL)
	{
		printf("No K_mer argument given. Please enter %s --help for more information\n", argv[0]);
		return 0;
	}
	
	int k = strlen(k_mer);
	char nucleicAcidKmer[k];
	char reverseComplementKmer[k];

	for (int i = 0; i < k; i++)
	{
		nucleicAcidKmer[i] = k_merLookupTable[k_mer[i]];
		reverseComplementKmer[k - i - 1] = reverseComplementTable[nucleicAcidKmer[i]];
		k_mer[i] = lookupTable[k_mer[i]];
	}

	char* dataArray = read_file(fileName, firstSequence, lookupTable);
	
	if (dataArray = NULL)
	{
		printf("Error reading from file: %s\n", fileName);
		exit(EXIT_FAILURE);
	}
	
	//print_sequences(firstSequence);
	
	currentSequence = firstSequence;

	if (chromosome == NULL)
	{
		while(currentSequence != NULL)
		{
			printf("checking %s\n", currentSequence->sequenceName);
			if (currentSequence->sequenceType == 1)
			{
				search_nucleic_acid(currentSequence, nucleicAcidKmer, reverseComplementKmer, k);
			} else if (currentSequence->sequenceType == 2)
			{
				search_amino_acid(currentSequence, k_mer, k);
			}	
			currentSequence = currentSequence->next;
		}
	} else
	{
		while(currentSequence != NULL)
		{
			if (strcmp(chromosome, currentSequence->sequenceName) == 0)
			{
				if (currentSequence->sequenceType == 1)
				{
					search_nucleic_acid(currentSequence, nucleicAcidKmer, reverseComplementKmer, k);
				} else if (currentSequence->sequenceType == 2)
				{
					search_amino_acid(currentSequence, k_mer, k);
				}
			}
			currentSequence = currentSequence->next;
		}
	}

	free(dataArray);
	
	return 0;
}


/*
*	Reads a file into memory and does preprocessing to aid in later searching.
*	Capable of loading files larger than 2GB.
*/
char* read_file(char* inputFile, Sequence* currentSequence, char* lookupTable)
{
	char* fileBuffer;
	int fileSize;
	int c;
	int read;
	int sequenceStart;
	int index = 0;
	
	struct stat64 statStruct;
	
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
	
	while((c = fgetc(file)) != EOF)
	{
		if (lookupTable[c] == 0)
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
				fileBuffer[index++] = 0;
			}
		} else
		{
			fileBuffer[index++] = lookupTable[c];
		}
	}
	
	currentSequence->sequenceLength = index - sequenceStart;

	return fileBuffer;
}

/*
*	Searches through a nucleic acid sequence and matchs the K-mer and the K-mer reverse complement.
*	Uses a binary AND operation to limit the number of conditionals within the for loop.
*/
void search_nucleic_acid(Sequence* currentSequence, char* k_mer, char* reverseComplementKmer, int k)
{
	int i;
	int j;
	float gcContent;
	int gc = 0;
	
	for (i = 0; i <= currentSequence->sequenceLength - k; i++)
	{
		if (k_mer[0] & currentSequence->data[i])
		{
			for (j = 1; j < k; j++)
			{
				if (!(k_mer[j] & currentSequence->data[ i + j ]))
				{
					break;
				}
			}
			if (j == k)
			{
				print_match(currentSequence, i, i + k, '+');
			}
		}
		
		if (reverseComplementKmer[0] & currentSequence->data[i])
		{
			for (j = 1; j < k; j++)
			{
				if (!(reverseComplementKmer[j] & currentSequence->data[ i + j ]))
				{
					break;
				}
			}
			if (j == k)
			{
				print_match(currentSequence, i, i + k, '-');
			}
		}
		if (currentSequence->data[i] & 6)
		{
			gc++;
		}
	}
	
	for (; i <= currentSequence->sequenceLength; i++)
	{
		if (currentSequence->data[i] & 6)
		{
			gc++;
		}
	}
	gcContent = 100 * ((float)gc / currentSequence->sequenceLength);
	printf("GC content for sequence %s: %.2f\%\n\n", currentSequence->sequenceName, gcContent);
}

/*
*	Searches through an amino acid sequence to match a given K-mer.
*	Uses exact matching to locate a match.
*/
void search_amino_acid(Sequence* currentSequence, char* k_mer, int k)
{
	int i;
	int j;
	
	for (i = 0; i < currentSequence->sequenceLength - k; i++)
	{
		if (k_mer[0] == currentSequence->data[i])
		{
			for (j = 1; j < k; j++)
			{
				if (k_mer[j] != currentSequence->data[i+j])
				{
					break;
				}
			}
			if (j == k)
			{
				print_match(currentSequence, i, i + k, '.');
			}
		}
	}
}

/*
*	Prints a found K-mer using GFF format.
*/
void print_match(Sequence* currentSequence, int start, int end, char strand)
{
	printf("%s \tMFW \tNo_Feature_Name \t%d \t%d \t. \t%c \t. No_GFF_Grouping_Attributes\n", currentSequence->sequenceName, start + 1, end, strand);
	fflush(stdout);
}


void print_sequences(Sequence* currentSequence)
{
	while (currentSequence != NULL)
	{
		printf("%s\n", currentSequence->sequenceName);
		fflush(stdout);
		
		for (int i = 0; i < currentSequence->sequenceLength; i++)
		{
			printf("%d, ", (int) currentSequence->data[i]);
		}
		printf("\n\n");
		
		currentSequence = currentSequence->next;
	}
}
