/*
* 	CSCI-5314 : HW4 : watson_hw4.c
*	Author: Mike Watson
*
*	Reads in a FASTA file and compares each sequence with each other sequence.
*	Builds a distance matrix and uses that to construct a MSA tree using Newick format
*
*	Sample command line: ./hw4 -f example.fasta -g 2 -s nuc.score
*
*/

#define __USE_LARGEFILE64
#define _LARGEFILE_SOURCE
#define _LARGEFILE64_SOURCE


#include <stdio.h>
#include <getopt.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include "watson_hw4.h"


int main (int argc, char** argv)
{
	char* fastaFile = NULL;
	char* scoringMatrixFile = NULL;
	int gapPenalty;
	char lookupTable[124] = {0};
	int numberOfSequences;
	
	Sequence* firstSequence = malloc(sizeof(Sequence));
	Sequence* currentSequence;
	
	/* -------Begin GETOPT parsing------- */
	
	int c;
	
	opterr = 0;
	
	while ((c = getopt (argc, argv, "f:g:s:")) != -1)
	{
		switch (c)
		{
			case 'f':
				/* handle fasta file argument */
				fastaFile = optarg;
				break;
			case 'g':
				/* handle global argument */
				gapPenalty = atoi(optarg);
			case 's':
				/* handle scoring matrix argument */
				scoringMatrixFile = optarg;
				break;
			case '?':
				if (optopt == 'f')
				{
					fastaFile = NULL;
				} else if (optopt == 's')
				{
					scoringMatrixFile = NULL;
				} else if (optopt == 'g')
				{
					gapPenalty = -19;
				} else
				{
					printf("Unknown option: %c\n", optopt);
				}
				return 1;
			default:
				abort ();
		}
	}
	
	/* -------End GETOPT parsing------- */
	
	if (fastaFile == NULL)
	{
		printf("Please specify a FASTA file for input. (use %s -f FILENAME to specify)\n", argv[0]);
		
		exit(1);
	}
	
	if (scoringMatrixFile == NULL)
	{
		printf("Please specify a Scoring Matrix file. (use %s -s FILENAME to specify)\n", argv[0]);
		
		exit(1);
	}
	
	if (gapPenalty < 0)
	{
		printf("Please specify a positive integer gap penalty. (use %s -g NUM to specify)\n", argv[0]);
	
		exit(1);
	}

	ScoringMatrix* score = read_scoring_file(scoringMatrixFile);
	
	char* fileBuffer = read_file(fastaFile, firstSequence, score, &numberOfSequences);
	
	if (fileBuffer == NULL)
	{
		/* error reading from file */
		printf ("Unable to read from file: %s. Exiting...", fastaFile);
		exit(0);
	}
	
	currentSequence = firstSequence;
	
	Node* root = build_tree(score, firstSequence, gapPenalty, numberOfSequences);

	print_tree(root);
	
	free(fileBuffer);
	
	return 0;
}


/*
*	Reads a scoring matrix from the specified file.
*
*	File should be of the following format:
*
*	A	3	-2	-1	-2
*	C	-2	3	-2	-1
*	G	-1	-2	3	-2
*	T	-2	-1	-2	3
*
*	Returns the scoring matrix as a ScoringMatrix struct containting the scoring matrix,
*	a lookup table for indexing into the matrix, and a reverse lookup table for converting
*	the numerical index back to a char.
*/
ScoringMatrix* read_scoring_file(char* scoreFile)
{
	ScoringMatrix* scoring = malloc(sizeof(ScoringMatrix));
	scoring->numChars = 0;
	
	char line[256];
	char* temp;
	int lineLength;
	int numLines = 0;
	int i;
	int j;
	
	FILE* file = fopen(scoreFile, "r");
	
	if (fgets(line, 256, file) != NULL)
	{
		strtok(line, " \t\n");
		
		while (strtok(NULL, " \t\n") != NULL)
		{
			scoring->numChars++;
		}
	} else
	{
		printf("Unable to read from score file. Exiting.\n");
		exit(1);
	}
	
	rewind(file);
	
	scoring->matrix = malloc(sizeof(int) * (scoring->numChars) * (scoring->numChars));
	scoring->lookupTable = calloc(128, sizeof(char));
	memset(scoring->lookupTable, -1, sizeof(char) * 128);

	scoring->lookupTableReverse = malloc (sizeof(char) * (scoring->numChars));
	
	for (i = 0; i < scoring->numChars; i++)
	{
		fgets(line, 256, file);

		temp = strtok(line, " \t\n");

		j = 0;
		scoring->lookupTable[temp[0]] = i;
		scoring->lookupTable[temp[0] ^ 32 ] = i;
		
		scoring->lookupTableReverse[i] = temp[0];
		
		while (((temp = strtok(NULL, " \t\n")) != NULL) && j <= scoring->numChars)
		{
			scoring->matrix[ i * (scoring->numChars) + j++] = atoi(temp);
		}
	}

	fclose(file);

	return scoring;
}


/*
*	Evaluates the distance between each sequence and builds a evolutionary tree.
*
*	ARGS
*		ScoringMatrix* score - ScoringMatrix struct used for scoring alignments
*		Sequence* currentSequence - the first sequence to start alignment from
*		int gapPenalty - the integer penalty for gaps
*		int numberOfSequences - number of sequences to include in the multiple sequnce alignment
*
*	RETURN
*		Node* - Node struct containing the top node of the resulting tree.
*
*/
Node* build_tree(ScoringMatrix* score, Sequence* currentSequence, int gapPenalty, int numberOfSequences)
{
	Sequence* otherSequence;
	Node* newNode;
	int maxDistance = 0;
	
	int i = 0;
	int j;
	
	float distanceMatrix[numberOfSequences][numberOfSequences];;
	Node* sequenceNodeArray[numberOfSequences];
	
	while (currentSequence != NULL)
	{
		j = i + 1;
		
		sequenceNodeArray[i] = malloc(sizeof(Node));
		sequenceNodeArray[i]->sequence = currentSequence;
		sequenceNodeArray[i]->height = 0;
		
		otherSequence = currentSequence->next;
		
		while (otherSequence != NULL)
		{
			distanceMatrix[i][i] = 0;
			distanceMatrix[i][j] = distanceMatrix[j][i] = global_alignment_distance(currentSequence, otherSequence, score, gapPenalty);
			
			if (distanceMatrix[i][j] > maxDistance)
			{
				maxDistance = distanceMatrix[i][j];
			}
			
			j++;
			otherSequence = otherSequence->next;
		}
		
		++i;
		currentSequence = currentSequence->next;
	}
	
	printf("Distance Matrix: (sequence names are omitted for ease of formatting)\n \t");
	for (int a = 0; a < numberOfSequences; a++)
	{
		printf("%c\t", 'a' + a);
	}
	printf("\n");
	
	for (int a = 0; a < numberOfSequences; a++)
	{
		printf("%c\t", 'a' + a);
		for (int b = 0; b < numberOfSequences; b++)
		{
			printf ("%.2f\t", distanceMatrix[a][b]);
		}
		printf("\n");
	}
	
	float minDistance;
	int min1;
	int min2;
	
	for (int a = 0; a < numberOfSequences - 1; a++)
	{
		minDistance = maxDistance + 1;
		
		for (i = 0; i < numberOfSequences; i++)
		{
			for (j = i + 1; j < numberOfSequences; j++)
			{
				if (distanceMatrix[i][j] > 0 && distanceMatrix[i][j] < minDistance)
				{
					minDistance = distanceMatrix[i][j];
					min1 = i;
					min2 = j;
				}
			}
		}
		
		sequenceNodeArray[min1]->distance = (minDistance / 2) - sequenceNodeArray[min1]->height;
		sequenceNodeArray[min2]->distance = (minDistance / 2) - sequenceNodeArray[min2]->height;
		
		newNode = malloc(sizeof(Node));
		newNode->sequence = NULL;
		newNode->distance = 0;
		newNode->left = sequenceNodeArray[min1];
		newNode->right = sequenceNodeArray[min2];
		newNode->height = sequenceNodeArray[min1]->distance + sequenceNodeArray[min1]->height;
		
		sequenceNodeArray[min1] = newNode;
		
		for (int i = 0; i < numberOfSequences; i++)
		{
			distanceMatrix[min1][i] = distanceMatrix[i][min1] = (distanceMatrix[min1][i] + distanceMatrix[i][min2]) / 2;
			distanceMatrix[i][min2] = distanceMatrix[min2][i] = -1;
		}
		distanceMatrix[min1][min1] = 0;
	}
	
	return newNode;
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
char* read_file(char* inputFile, Sequence* currentSequence, ScoringMatrix* score, int* numberOfSequences)
{
	char* fileBuffer;
	int fileSize;
	int c;
	int read;
	int sequenceStart;
	int index = 0;
	*numberOfSequences = 0;
	
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
			*numberOfSequences += 1;
			currentSequence->sequenceName = &fileBuffer[index];
		} else if (c == ' ' || c == '\t')
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
		if (score->lookupTable[c] < 0)
		{
			if (c == '>')
			{
				*numberOfSequences += 1;
				currentSequence->sequenceLength = index - sequenceStart;
				currentSequence->next = malloc(sizeof(Sequence));
				currentSequence = currentSequence->next;
				currentSequence->sequenceName = &fileBuffer[index];
				currentSequence->sequenceType = 1;
				while ( (c = fgetc(file)) != '\n')
				{
					if (c == ' ' || c == '\t')
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
			} else if (c == '\n')
			{
			}else
			{
				printf("Error: character '%c' is not defined in the scoring matrix file. Exiting.\n", c);
				exit(1);
			}
		} else
		{
			fileBuffer[index++] = score->lookupTable[c];
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
*		ScoringMatrix* score - ScoringMatrix struct containing the scoring information for evaluating matches.
*		int gapPenalty - the integer penalty applied for any gap.
*
*	RETURN:  float: distance score between the two input sequences. Evaluated as a distance between 0 - 50.
*					Any negative alignment scores will result in a distance higher than 50.
*/
float global_alignment_distance(Sequence* sequenceOne, Sequence* sequenceTwo, ScoringMatrix* score, int gapPenalty)
{
	AlignmentMatrix* alignment = malloc(sizeof(AlignmentMatrix));
	
	alignment->matrixWidth = sequenceOne->sequenceLength + 1;
	alignment->matrixHeight = sequenceTwo->sequenceLength + 1;
	int bestScore = 0;
	
	int directionTable[3] = { 1, alignment->matrixWidth + 1, alignment->matrixWidth };
	
	alignment->score = malloc(alignment->matrixWidth * alignment->matrixHeight * sizeof(int));

	int i;
	int j;
	int index;
	int left;
	int diag;
	int up;
	
	alignment->score[0] = 0;
	
	for (j = 1; j < alignment->matrixWidth; j++)
	{
		alignment->score[j] = j * GLOBAL_GAP;
		bestScore += score->matrix[sequenceOne->data[j-1] * (score->numChars) + sequenceOne->data[j-1]];
	}
	
	for (i = 1; i < alignment->matrixHeight; i++)
	{
		alignment->score[i * alignment->matrixWidth] = i * -gapPenalty;
		for (j = 1; j < alignment->matrixWidth; j++)
		{
			index = i * alignment->matrixWidth + j;
			left = alignment->score[index - directionTable[LEFT]] - gapPenalty;
			diag = alignment->score[index - directionTable[DIAG]]
								+ score->matrix[sequenceOne->data[j-1] * (score->numChars) + sequenceTwo->data[i-1]];
			up = alignment->score[index - directionTable[UP]] - gapPenalty;
			
			if ( left > diag )
			{
				if ( left > up )
				{
					alignment->score[index] = left;
				}
				else
				{
					alignment->score[index] = up;
				}
			}
			else
			{
				if ( up > diag )
				{
					alignment->score[index] = up;
				}
				else
				{
					alignment->score[index] = diag;
				}
			}
		}
	}
	
	
	float returnVal = ((bestScore - (float)alignment->score[index]) * 50)/bestScore;

	free(alignment->score);
	free(alignment);
	
	return returnVal ;
}

/*
*	Uses a recursive function to output the tree result in Newick format.
*
*/
void print_tree(Node* root)
{
	printf("\n");
	print_tree_recursive(root);
	
	printf(";\n");
}


/*
*	Recursive function used in printing the Newick format tree.
*/
void print_tree_recursive(Node* aNode)
{
	if (aNode->sequence != NULL)
	{
		printf("%s:%.2f", aNode->sequence->sequenceName, aNode->distance);
	} else
	{
		printf("(");
		print_tree_recursive(aNode->left);
		printf(",");
		print_tree_recursive(aNode->right);
		printf("):%.2f", aNode->distance);
	}
}