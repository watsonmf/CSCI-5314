/*
*
*
*
*
*/

int main (int argc, char** argv)
{
	/*	-------Begin POPT parsing-------	*/

	poptContext	POPT_Context;  /* context for parsing command-line options */
	int			POPT_Ret;      /* used for iterating over the arguments */

	struct poptOption optionsTable[] =
	{
		{ "file", 'f', POPT_ARG_STRING, &fileName, 'f', "Specify FASTA format genome file to read in", "FILENAME" },
		{ "chromosome", 'c', POPT_ARG_STRING, &chromosome, 'c', "Specify chromosome to search (default is all)", "NAME" },
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

		}
	}
	
	fflush(stdout);

	if (POPT_Ret < -1)
	{
		/* an error occurred during option processing */
		fprintf(stderr, "%s: %s\n",
		        poptBadOption(POPT_Context, POPT_BADOPTION_NOALIAS),
		        poptStrerror(POPT_Ret));
		return 1;
	}
	/*	-------End POPT parsing-------	*/
}



