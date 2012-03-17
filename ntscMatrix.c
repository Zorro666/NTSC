#include <malloc.h>
#include <stdio.h>

float** matrixMalloc(const unsigned int numRows, const unsigned int numCols)
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
	return result;
}

void matrixFree(float** matrix, const unsigned int numRows)
{
	unsigned int i;
	for (i = 0; i < numRows; i++)
	{
		free(matrix[i]);
	}
	free(matrix);
}

#if MATRIX_DEBUG
void matrixPrintf(float** mat, unsigned int numRows, unsigned int numCols, const char* const name)
{
	unsigned int i;
	for (i = 0; i < numRows; i++)
	{
		unsigned int j;
		for (j = 0; j < numCols; j++)
		{
			printf("%s[%d][%d] = %f ", name, i, j, mat[i][j]);
		}
		printf("\n");
	}
}
#endif /* #if MATRIX_DEBUG */

void matrixMultiply(float** result, float** left, float** right, 
									  const unsigned int numRowsLeft, const unsigned int numColsLeft, 
										const unsigned int numColsRight)
{
	unsigned int i;
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

void matrixTranspose(float** transpose, float** matrix, const unsigned int numRows, const unsigned int numCols)
{
	unsigned int i;
	for (i = 0; i < numRows; i++)
	{
		unsigned int j;
		for (j = 0; j < numCols; j++)
		{
			transpose[j][i] = matrix[i][j];
		}
	}
}

