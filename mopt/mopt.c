


typedef struct _OPTION
{
	char* longName;
	char shortName;
	int argInfo;
	void * arg;
	int val;
	char* description;
	char* argDescription;
} moptOption;




MOPT_Context moptGetContext(NULL, int argc, const char** argv, moptOption, int question)
{
		
}

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