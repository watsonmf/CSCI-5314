/*
* 	CSCI-5314 : HW3 : watson_hw3.c
*	Author: Mike Watson
*
*	Reads in a FASTA file and compares each pair of sequences to find either the global or local alignment.
*
*	This program uses a hash table to store the number of unique k-mers found in all the sequences.
*
*	Required libraries: POPT
*
*	Sample command line: ./hw3 -f example.fasta -g
*
*
*	In the event that the program is run with the -l option for local alignment, it will check each
*	Sequence pair to find the shortest sequence and run the alignment with that as the bottom sequence.
*	Trailing gaps are not penalized for the shorter sequence (bottom sequence).
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

#include "watson_hw3.h"


int main (int argc, char** argv)
{
	char* fileName = NULL;
	int globalAlignment = 1;
	int localAlignment = 0;
	
	Sequence* firstSequence = malloc(sizeof(Sequence));
	Sequence* currentSequence;
	
	/*	-------Begin POPT parsing-------	*/

	poptContext	POPT_Context;  /* context for parsing command-line options */
	int			POPT_Ret;      /* used for iterating over the arguments */

	struct poptOption optionsTable[] =
	{
		{ "file", 'f', POPT_ARG_STRING, &fileName, 'f', "Specify FASTA format genome file to read in", "FILENAME" },
		{ "global", 'g', POPT_ARG_NONE, NULL, 'g', "Use global alignment algorithm", 0 },
		{ "local", 'l', POPT_ARG_NONE, NULL, 'l', "Use local alignment algorithm", 0 },
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
		case 'g':
			/* handle global alignment argument */
			globalAlignment = 1;
			localAlignment = 0;
			break;
		case 'l':
			localAlignment = 1;
			globalAlignment = 0;
			/* handle local alignment argument */
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
	
	if (globalAlignment == 0 && localAlignment == 0)
	{
		printf("Please specify either global or local alignment algorithms. (%s --help to list all options)\n", argv[0]);
	
		exit(1);
	}

	char* fileBuffer = read_file(fileName, firstSequence);
	
	if (fileBuffer == NULL)
	{
		/* error reading from file */
		printf ("Unable to read from file: %s. Exiting...", fileName);
		exit(0);
	}
	
	currentSequence = firstSequence;
	AlignmentMatrix* aligned;
	
	if (globalAlignment) // Use global alignment algorithm
	{
		printf("Calculating global alignment for %s:\n\n", fileName);
		
		// iterate through the sequence linked list two at a time and compair each pair of sequences.
		while (currentSequence != NULL)
		{
			if (currentSequence->next == NULL)
			{
				printf("File contained an odd number of sequences, last sequence will not be compared.\n");
				break;
			}
			
			AlignmentMatrix* aligned = global_alignment(currentSequence, currentSequence->next);
			print_global_alignment(aligned, currentSequence, currentSequence->next);
			
			currentSequence = currentSequence->next->next;
			free(aligned);
		}
	}
	
	if (localAlignment)// Use local alignment algorithm
	{
		printf("Calculating local alignment for %s:\n\n", fileName);
		
		// iterate through the sequence linked list two at a time and compair each pair of sequences.
		while (currentSequence != NULL)
		{
			if (currentSequence->next == NULL)
			{
				printf("File contained an odd number of sequences, last sequence will not be compared.\n");
				break;
			}
			
			// always match the shorter sequence on the bottom
			if (currentSequence->sequenceLength >= currentSequence->next->sequenceLength)
			{
				aligned = local_alignment(currentSequence, currentSequence->next);
			} else
			{
				aligned = local_alignment(currentSequence->next, currentSequence);
			}
			
			print_local_alignment(aligned, currentSequence, currentSequence->next);
			
			currentSequence = currentSequence->next->next;
			free(aligned);
		}
	}
	
	free(fileBuffer);
	
	return 0;
}

/*
*	Reads a file into memory and builds a linked-list containing all the sequences found in the file.
*	Uses stat64() and fopen64() to allow for loading files larger than 2GB.
*
*	All file data is stored in the fileBuffer array and each of the sequence structs use pointers to
*	access their respective data.
*
*	ARGS
*		char* inputFile - String containing the file location.
*		Sequence* currentSequence - pointer to the memory location of the first member of a linked list
*									that will contain all the sequences read in from the file.
*
*	RETURN: char* - returns a pointer to an array containing all data read in from the file. The sequence
					linked list contains pointers to their respective information.
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
*	Calculates the global alignment using the Needleman-Wunsch algorithm
*
*	ARGS:
*		Sequence* sequenceOne - horizonatal sequence used for alignemnt
*		Sequence* sequenceTwo - verticle sequence used for alignment
*
*	RETURN:  AlignmentMatrix* - Returns an AlignmentMatrix struct containing the 2d dynamic programming
*					matrices for score and direction.
*/
AlignmentMatrix* global_alignment(Sequence* sequenceOne, Sequence* sequenceTwo)
{
	AlignmentMatrix* alignment = malloc(sizeof(AlignmentMatrix));
	
	alignment->matrixWidth = sequenceOne->sequenceLength + 1;
	alignment->matrixHeight = sequenceTwo->sequenceLength + 1;
	
	int directionTable[3] = { 1, alignment->matrixWidth + 1, alignment->matrixWidth };
	
	alignment->score = malloc(alignment->matrixWidth * alignment->matrixHeight * sizeof(int));
	alignment->direction = malloc(alignment->matrixWidth * alignment->matrixHeight * sizeof(char));
	
	int i;
	int j;
	int index;
	int left;
	int diag;
	int up;
	
	for (j = 0; j < alignment->matrixWidth; j++)
	{
		alignment->score[j] = j * GLOBAL_GAP;
		alignment->direction[j] = LEFT;
	}
	
	for (i = 1; i < alignment->matrixHeight; i++)
	{
		alignment->score[i * alignment->matrixWidth] = i * GLOBAL_GAP;
		alignment->direction[i * alignment->matrixWidth] = UP;
		
		for (j = 1; j < alignment->matrixWidth; j++)
		{
			index = i * alignment->matrixWidth + j;
			left = alignment->score[index - directionTable[LEFT]] + GLOBAL_GAP;
			diag = alignment->score[index - directionTable[DIAG]]
								+	(sequenceOne->data[j - 1] == sequenceTwo->data[i - 1] ? GLOBAL_MATCH : GLOBAL_MISMATCH);
			up = alignment->score[index - directionTable[UP]] + GLOBAL_GAP;
			
			if ( left < diag )
			{
				if ( left < up )
				{
					alignment->score[index] = left;
					alignment->direction[index] = LEFT;
				}
				else
				{
					alignment->score[index] = up;
					alignment->direction[index] = UP;
				}
			}
			else
			{
				if ( up < diag )
				{
					alignment->score[index] = up;
					alignment->direction[index] = UP;
				}
				else
				{
					alignment->score[index] = diag;
					alignment->direction[index] = DIAG;
				}
			}
		}
	}
	
	return alignment;
}

/*
*	Calculates the local alignment using the Smith-Waterman algorithm.
*
*	WARNING: It is assumed that sequenceTwo is the shorter sequence than sequenceOne. Trailing gaps
*			on sequenceTwo are free.
*
*	ARGS:
*		Sequence* sequenceOne - horizonatal sequence used for alignemnt
*		Sequence* sequenceTwo - verticle sequence used for alignment
*
*	RETURN:  AlignmentMatrix* - Returns an AlignmentMatrix struct containing the 2d dynamic programming
*					matrices for score and direction.
*/
AlignmentMatrix* local_alignment( Sequence* sequenceOne, Sequence* sequenceTwo )
{
	AlignmentMatrix* alignment = malloc(sizeof(AlignmentMatrix));
	
	alignment->matrixWidth = sequenceOne->sequenceLength + 1;
	alignment->matrixHeight = sequenceTwo->sequenceLength + 1;
	
	int directionTable[3] = { 1, alignment->matrixWidth + 1, alignment->matrixWidth };
	
	alignment->score = malloc(alignment->matrixWidth * alignment->matrixHeight * sizeof(int));
	alignment->direction = malloc(alignment->matrixWidth * alignment->matrixHeight * sizeof(char));
	
	int i;
	int j;
	int index;
	int left;
	int diag;
	int up;
	
	// Initialize the first row to zero to represent no charge for leading gap
	for (j = 0; j < alignment->matrixWidth; j++)
	{
		alignment->score[j] = 0;
		alignment->direction[j] = LEFT;
	}
	
	// calculate the score the rest of the matrix except the last row
	for (i = 1; i < alignment->matrixHeight - 1; i++)
	{
		alignment->score[i * alignment->matrixWidth] = 0;
		alignment->direction[i * alignment->matrixWidth] = UP;
		
		for (j = 1; j < alignment->matrixWidth; j++)
		{
			index = i * alignment->matrixWidth + j;
			left = alignment->score[index - directionTable[LEFT]] + LOCAL_GAP;
			diag = alignment->score[index - directionTable[DIAG]]
								+	(sequenceOne->data[j - 1] == sequenceTwo->data[i - 1] ? LOCAL_MATCH : LOCAL_MISMATCH);
			up = alignment->score[index - directionTable[UP]] + LOCAL_GAP;
			
			if ( left > diag )
			{
				if ( left > up )
				{
					alignment->score[index] = left;
					alignment->direction[index] = LEFT;
				}
				else
				{
					alignment->score[index] = up;
					alignment->direction[index] = UP;
				}
			}
			else
			{
				if ( up > diag )
				{
					alignment->score[index] = up;
					alignment->direction[index] = UP;
				}
				else
				{
					alignment->score[index] = diag;
					alignment->direction[index] = DIAG;
				}
			}
			
			if (alignment->score[index] < 0)
			{
				alignment->score[index] = 0;
			}
		}
	}
	
	index++;
	
	// calculate the score for the last row to allow for no charge for trailing gaps on the bottom sequence.
	for (j = 1; j < alignment->matrixWidth; j++)
	{
		index++;
		left = alignment->score[index - directionTable[LEFT]] + 0;
		diag = alignment->score[index - directionTable[DIAG]]
							+	(sequenceOne->data[j - 1] == sequenceTwo->data[i - 1] ? LOCAL_MATCH : LOCAL_MISMATCH);
		up = alignment->score[index - directionTable[UP]] + LOCAL_GAP;
		
		if ( left > diag )
		{
			if ( left > up )
			{
				alignment->score[index] = left;
				alignment->direction[index] = LEFT;
			}
			else
			{
				alignment->score[index] = up;
				alignment->direction[index] = UP;
			}
		}
		else
		{
			if ( up > diag )
			{
				alignment->score[index] = up;
				alignment->direction[index] = UP;
			}
			else
			{
				alignment->score[index] = diag;
				alignment->direction[index] = DIAG;
			}
		}
	}

	return alignment;
}

/*
*	Prints the results from the global_alignment function. Backtracks through the 2d alignment matrix,
*		starting at the last index backtracking through until the index is 0. Builds a string for the
*		top and bottom aligned sequences and then prints the result to stdout.
*
*	ARGS:
*		AlignmentMatrix* alignment - AlignmentMatrix struct containing the 2d dynamic programming
*					matrices for score and direction.
*		Sequence* sequenceOne - horizonatal sequence used for alignemnt
*		Sequence* sequenceTwo - verticle sequence used for alignment
*/
void print_global_alignment(AlignmentMatrix* alignment, Sequence* sequenceOne, Sequence* sequenceTwo)
{
	int index = alignment->matrixWidth * alignment->matrixHeight - 1;
	int score = alignment->score[index];
	
	// Set alignedSeuqenceIndex to be the max possible length of the aligned sequence (plus one for string end marker '\0'.
	int alignedSequenceIndex = alignment->matrixWidth + alignment->score[alignment->matrixWidth * alignment->matrixHeight - 1] + 1;
	
	int sequenceOneIndex = sequenceOne->sequenceLength - 1;
	int sequenceTwoIndex = sequenceTwo->sequenceLength - 1;
	
	int directionTable[3] = { 1, alignment->matrixWidth + 1, alignment->matrixWidth };
	
	char alignedSequenceOne[alignedSequenceIndex];
	char alignedSequenceTwo[alignedSequenceIndex];
	
	alignedSequenceIndex--;
	
	alignedSequenceOne[alignedSequenceIndex] = '\0';
	alignedSequenceTwo[alignedSequenceIndex--] = '\0';
	
	while (index > 0)
	{
		switch ( alignment->direction[index] )
		{
		case LEFT:
			alignedSequenceOne[alignedSequenceIndex] = sequenceOne->data[sequenceOneIndex--];
			alignedSequenceTwo[alignedSequenceIndex] = '-';
			break;
		case DIAG:
			alignedSequenceOne[alignedSequenceIndex] = sequenceOne->data[sequenceOneIndex--];
			alignedSequenceTwo[alignedSequenceIndex] = sequenceTwo->data[sequenceTwoIndex--];
			break;
		case UP:
			alignedSequenceOne[alignedSequenceIndex] = '-';
			alignedSequenceTwo[alignedSequenceIndex] = sequenceTwo->data[sequenceTwoIndex--];
			break;
		}
		
		alignedSequenceIndex--;
		index -= directionTable[alignment->direction[index]];
	}
	
	alignedSequenceIndex++;
	
	printf("%s\n%s\n", &alignedSequenceOne[alignedSequenceIndex], &alignedSequenceTwo[alignedSequenceIndex]);
	
	for (int i = 0; i < strlen(&alignedSequenceTwo[alignedSequenceIndex]); i++)
	{
		printf("=");
	}
	
	printf(" (%d)\n", score);
}


/*
*	Prints the results from the local_alignment function. Backtracks through the 2d alignment matrix,
*		starting at the last index backtracking through until the index is 0. Builds a string for the
*		top and bottom aligned sequences and then prints the result to stdout.
*
*	ARGS:
*		AlignmentMatrix* alignment - AlignmentMatrix struct containing the 2d dynamic programming
*					matrices for score and direction.
*		Sequence* sequenceOne - horizonatal sequence used for alignemnt
*		Sequence* sequenceTwo - verticle sequence used for alignment
*/
void print_local_alignment(AlignmentMatrix* alignment, Sequence* sequenceOne, Sequence* sequenceTwo)
{
	int index = alignment->matrixWidth * alignment->matrixHeight - 1;
	int score = 0;
	
	int sequenceOneIndex = sequenceOne->sequenceLength - 1;
	int sequenceTwoIndex = sequenceTwo->sequenceLength - 1;
	
	int directionTable[3] = { 1, alignment->matrixWidth + 1, alignment->matrixWidth };
	
	char alignedSequenceOne[sequenceOne->sequenceLength + sequenceTwo->sequenceLength];
	char alignedSequenceTwo[sequenceOne->sequenceLength + sequenceTwo->sequenceLength];
	
	int alignedSequenceIndex = sequenceOne->sequenceLength + sequenceTwo->sequenceLength - 1;
	
	alignedSequenceOne[alignedSequenceIndex] = '\0';
	alignedSequenceTwo[alignedSequenceIndex--] = '\0';
	
	while (sequenceOneIndex >=0 || sequenceTwoIndex >= 0)
	{
		switch ( alignment->direction[index] )
		{
		case LEFT:
			// trailing gaps are not counted on sequence two
			if (sequenceTwoIndex >= 0 && sequenceTwoIndex < sequenceTwo->sequenceLength - 1)
			{
				score++;
			}
			alignedSequenceOne[alignedSequenceIndex] = sequenceOne->data[sequenceOneIndex--];
			alignedSequenceTwo[alignedSequenceIndex] = '-';
			break;
		case DIAG:
			if (sequenceOne->data[sequenceOneIndex] != sequenceTwo->data[sequenceTwoIndex])
			{
				score++;
			}
			alignedSequenceOne[alignedSequenceIndex] = sequenceOne->data[sequenceOneIndex--];
			alignedSequenceTwo[alignedSequenceIndex] = sequenceTwo->data[sequenceTwoIndex--];
			break;
		case UP:
			if (sequenceOneIndex >= 0)
			{
				score++;
			}
			alignedSequenceOne[alignedSequenceIndex] = '-';
			alignedSequenceTwo[alignedSequenceIndex] = sequenceTwo->data[sequenceTwoIndex--];
			break;
		}
		
		alignedSequenceIndex--;
		index -= directionTable[alignment->direction[index]];
	}
	
	alignedSequenceIndex++;
	
	printf("%s\n%s\n", &alignedSequenceOne[alignedSequenceIndex], &alignedSequenceTwo[alignedSequenceIndex]);
	
	for (int i = 0; i < strlen(&alignedSequenceTwo[alignedSequenceIndex]); i++)
	{
		printf("=");
	}
	
	printf(" (%d)\n", score);
}

/*
* Test function used for printing the dynamic programming matrix
*/
void print_matrix(AlignmentMatrix* alignment)
{
	char directionTable[3] = {'l', 'd', 'u'};
	int index;
	
	for (int i = 0; i < alignment->matrixHeight; i++)
	{
		for (int j = 0; j < alignment->matrixWidth; j++)
		{
			index = i * alignment->matrixWidth + j;
			printf("%c%d ", directionTable[alignment->direction[index]], alignment->score[index]);
		}
		printf("\n");
	}
	printf("\n");
}




