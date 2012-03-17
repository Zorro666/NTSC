#include <malloc.h>
#include <stdio.h>
#include <stdlib.h>

#include "ntscMatrix.h"

void matrixCreate(NtscMatrix* const ntscMatrix, const unsigned int numRows, const unsigned int numCols, const char* const name)
{
	unsigned int i;
	float** result = NULL;
	result = malloc(sizeof(float*)*numRows);
	for (i = 0; i < numRows; i++)
	{
		unsigned int j;
		result[i] = malloc(sizeof(float)*numCols);
		for (j = 0; j < numCols; j++)
		{
			result[i][j] = 0.0f;
		}
	}
	ntscMatrix->m_matrix = result;
	ntscMatrix->m_numRows = numRows;
	ntscMatrix->m_numCols = numCols;
	ntscMatrix->m_name = name;
}

void matrixFree(NtscMatrix* const ntscMatrix)
{
	unsigned int i;
	const unsigned int numRows = ntscMatrix->m_numRows;
	float** const matrix = ntscMatrix->m_matrix;
	for (i = 0; i < numRows; i++)
	{
		free(matrix[i]);
	}
	free(matrix);

	ntscMatrix->m_matrix = NULL;
	ntscMatrix->m_numRows = 0;
	ntscMatrix->m_numCols = 0;
}

#if MATRIX_DEBUG
void matrixPrintf(const NtscMatrix* const ntscMatrix, const char* const inputName)
{
	unsigned int i;
	const unsigned int numRows = ntscMatrix->m_numRows;
	const unsigned int numCols = ntscMatrix->m_numCols;
	const float** const matrix = (const float** const)ntscMatrix->m_matrix;
	const char* const name = (inputName != NULL) ? inputName : ntscMatrix->m_name;
	for (i = 0; i < numRows; i++)
	{
		unsigned int j;
		for (j = 0; j < numCols; j++)
		{
			printf("%s[%d][%d] = %f ", name, i, j, matrix[i][j]);
		}
		printf("\n");
	}
}
#endif /* #if MATRIX_DEBUG */

/* result = left * right : left = m x n, right = n x o, result = m x o */
void matrixMultiply(NtscMatrix* const ntscResult, const NtscMatrix* const ntscLeft, const NtscMatrix* const ntscRight)
{
	unsigned int i;
	const unsigned int numRowsLeft = ntscLeft->m_numRows;
	const unsigned int numColsLeft = ntscLeft->m_numCols;
	const unsigned int numColsRight = ntscRight->m_numCols;
	const float** const left = (const float** const)ntscLeft->m_matrix;
	const float** const right = (const float** const)ntscRight->m_matrix;
	float** const result = ntscResult->m_matrix;

	if (numRowsLeft != ntscResult->m_numRows)
	{
		fprintf(stderr,"ERROR: matrixMultiply: left:'%s' result:'%s' numRowsLeft != numRowsResult %d != %d\n",
						ntscLeft->m_name, ntscResult->m_name, numRowsLeft, ntscResult->m_numRows);
		exit(-1);
	}
	if (numColsRight != ntscResult->m_numCols)
	{
		fprintf(stderr,"ERROR: matrixMultiply: right:'%s' result:'%s' numColsRight != numColsResult %d != %d\n",
						ntscRight->m_name, ntscResult->m_name, numColsRight, ntscResult->m_numCols);
		exit(-1);
	}
	if (numColsLeft != ntscRight->m_numRows)
	{
		fprintf(stderr,"ERROR: matrixMultiply: left:'%s' right:'%s' numColsLeft != numRowsRight %d != %d\n",
						ntscLeft->m_name, ntscRight->m_name, numRowsLeft, ntscRight->m_numRows);
		exit(-1);
	}
	for (i = 0; i < numRowsLeft; i++)
	{
		unsigned int j;
		for (j = 0; j < numColsRight; j++)
		{
			unsigned int k;
			result[i][j] = 0.0f;
			for (k = 0; k < numColsLeft; k++)
			{
				result[i][j] += left[i][k] * right[k][j];
			}
		}
	}
}

void matrixTranspose(NtscMatrix* const ntscTranspose, const NtscMatrix* const ntscMatrix)
{
	unsigned int i;
	const unsigned int numRows = ntscMatrix->m_numRows;
	const unsigned int numCols = ntscMatrix->m_numCols;
	const float** const matrix = (const float** const)ntscMatrix->m_matrix;
	float** const transpose = ntscTranspose->m_matrix;
	if (ntscTranspose->m_numRows != numRows)
	{
		fprintf(stderr,"ERROR: matrixTranpose: transpose:'%s' matrix:'%s' numRows doesn't match %d != %d\n",
						ntscTranspose->m_name, ntscMatrix->m_name, ntscTranspose->m_numRows, numRows);
		exit(-1);
	}
	if (ntscTranspose->m_numCols != numCols)
	{
		fprintf(stderr,"ERROR: matrixTranpose: transpose:'%s' matrix:'%s' numCols doesn't match %d != %d\n",
						ntscTranspose->m_name, ntscMatrix->m_name, ntscTranspose->m_numCols, numCols);
		exit(-1);
	}
	for (i = 0; i < numRows; i++)
	{
		unsigned int j;
		for (j = 0; j < numCols; j++)
		{
			transpose[j][i] = matrix[i][j];
		}
	}
}

