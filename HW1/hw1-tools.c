
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>



/*
*	
*
*/
char* read_file(char* inputFile, sequence* currentSequence)
{
	char* fileBuffer;
	int fileSize;
	int c;
	int read;
	int total = 0;
	
	struct stat statStruct;
	
	if (stat64(inputFile, &statStruct) == -1) 
	{
		return NULL;
	}
	
	fileSize = statStruct.st_size;
	
	FILE* file = fopen64(inputFile, "r");

	fileBuffer = malloc(fileSize);
	
	while((c = fgetc(file)) != EOF)
	{
		if (lookupTable[c] == 0)
		{
			if (c == '>')
			{
				currentSequence->next = malloc(sizeof(sequence));
				currentSequence = currentSequence->next;
				currentSequence->fileName = fileBuffer[++i];
				while ( (c = fgetc(file)) != '\n')
				{
					if (c == ' ')
					{
						fileBuffer[index++] = '\0';
						currentSequence->sequenceDescription = fileBuffer[index];
					} else
					{
						fileBuffer[index++] = (char) c;
					}
				}
				
				fileBuffer[index++] = '\0';
				currentSequence->data = fileBuffer[index];
				
			} 
		} else
		{
			fileBuffer[index++] = lookupTable[c];
		}
	}

	return fileBuffer;
}

void search_nucleic_acid(sequence* currentSequence, char* k_mer, char* reverseComplementKmer)
{
	int i;
	int j;
	int k = strlen(k_mer);
	char reverseComplementKmer[k];
	
	for (i = 0; i < k; i++)
	{
		reverseComplementKmer[i] = reverseComplementTable[k_mer[k - i - 1]];
	}
	
	for (i = 0; i < currentSequence->sequenceLength - k; i++)
	{
		if (k_mer[0] & currentSequence->data[i])
		{
			for (j = 1; j < k; j++)
			{
				if (!(k_mer[j] & currentSequence->data[i+j]))
				{
					break;
				}
			}
			if (j == k)
			{
				print_match(currentSequence, i, '+');
			}
		}
	}
}

void search_amino_acid(sequence* currentSequence, char* k_mer)
{
	
}


void print_match(sequence* currentSequence, int location, char strand)
{
	
}