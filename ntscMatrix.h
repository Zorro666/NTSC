#ifndef NTSCMATRIX_H
#define NTSCMATRIX_H

#define MATRIX_DEBUG 0

typedef struct NtscMatrix
{
	float** m_matrix;
	unsigned int m_numRows;
	unsigned int m_numCols;
	const char* m_name;
} NtscMatrix;

void matrixCreate(NtscMatrix* const ntscMatrix, const unsigned int numRows, const unsigned int numCols, const char* const name);
void matrixFree(NtscMatrix* const ntscMatrix);
void matrixPrintf(const NtscMatrix* const ntscMatrix, const char* const inputName);
void matrixMultiply(NtscMatrix* const ntscResult, const NtscMatrix* const ntscLeft, const NtscMatrix* const ntscRight);
void matrixTranspose(NtscMatrix* const ntscTranspose, const NtscMatrix* const ntscMatrix);

#endif /* #ifndef NTSCMATRIX_H */

