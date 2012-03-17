#ifndef NTSCMATRIX_H
#define NTSCMATRIX_H

#define MATRIX_DEBUG 0

typedef struct ntscMatrix
{
	float** m_matrix;
	int m_numRows;
	int m_numCols;
} ntscMatrix;

float** matrixMalloc(const unsigned int numRows, const unsigned int numCols);

void matrixFree(float** matrix, const unsigned int numRows);

void matrixPrintf(float** mat, unsigned int numRows, unsigned int numCols, const char* const name);

void matrixMultiply(float** result, float** left, float** right, 
									  const unsigned int numRowsLeft, const unsigned int numColsLeft, 
										const unsigned int numColsRight);

void matrixTranspose(float** transpose, float** matrix, const unsigned int numRows, const unsigned int numCols);

#endif /* #ifndef NTSCMATRIX_H */

