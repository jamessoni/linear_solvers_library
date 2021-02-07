#pragma once
#include "Matrix.h"

template <class T>
class CSRMatrix : public Matrix<T>
{
public:

	// constructor where we want to preallocate ourselves
	CSRMatrix(int rows, int cols, int nnzs, bool preallocate);
	// constructor where we already have allocated memory outside
	CSRMatrix(int rows, int cols, int nnzs, T* values_ptr, int* row_position, int* col_index);
	// destructor
	~CSRMatrix();

	// Print out the values in our matrix
	virtual void printMatrix();

	// Perform some operations with our matrix
	void matVecMult(T* input, T* output);
	void dense2sparse(Matrix<T>& tosparsify, CSRMatrix<T>* output);

	void jacobi_solver_sparse(CSRMatrix<T>* A, T* b, T* x, int maxIter, bool initialised, float tol);

	float RMS_norm_diff(T* vec_a, T* vec_b);

	void gauss_seidel_sparse(CSRMatrix<T>* A, T* b, T* x, int maxIter, float tol);

  	CSRMatrix<T>* matMatMult(CSRMatrix<T>& mat_right);
	
	void transposeiflower();
	CSRMatrix<T>* CholeskyDecomp();
	void fsubstitution(T* b, T* y);
	void bsubstitution(T* y, T* x);
	void CholeskySolve(Matrix<T>* A, T* b, T* x);


	// Explicitly using the C++11 nullptr here
	int* row_position = nullptr;
	int* col_index = nullptr;

	// How many non-zero entries we have in the matrix
	int nnzs = -1;

	// Private variables - there is no need for other classes 
	// to know about these variables
private:
};