/*
*
*	A 00000001
*	C 00000010
*	G 00000100
*	T 00001000
*	U 00010000 U not needed
*	
*
*
*
*
*
*
*
*
*
*
*
*
*
*/


#define LOOKUP_TABLE {0,0,0,0,0,0,0,0,0,0,'\n',0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,   \
!,",#,$, %, &, ', (, ), *, +, ,, -, ., /, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, :, ;, <, =, >, ?, @, \
A, B, C, D, E, F, G, H, I, J, K, L, M, N, O, P, Q, R, S, T, U, V, W, X, Y, Z, [, \, ], ^, _, `,  \
a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p, q, r, s, t, u, v, w, x, y, z} 

typedef char type;

typedef struct _SEQUENCE
{
	struct _SEQUENCE* next;
	char* sequenceName;
	char* sequenceDescription;
	int sequenceLength;
	type sequenceType; // 1 = nucleic acid, 2 = amino acid
	char* data; // pointer to the start of sequence.
} sequence;

typedef struct _FOUNDLIST // linked list for matched patterns
{
	
};


