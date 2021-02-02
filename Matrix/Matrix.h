#pragma once

template <class T>
class Matrix
{
public:

   // constructor where we want to preallocate ourselves
   Matrix(int rows, int cols, bool preallocate);
   // constructor where we already have allocated memory outside
   Matrix(int rows, int cols, T *values_ptr);
   // destructor
   virtual ~Matrix();

   // Print out the values in our matrix
   void printValues();
	virtual void printMatrix();

   // Perform some operations with our matrix
   void matMatMult(Matrix<T>& mat_left, Matrix<T>& output);
   void matVecMult(T* vec, T* output);
   void vecVecsubtract(T* vec_a, T* vec_b, T* output);
   float RMS_norm_diff(T* vec_a, T* vec_b);

   // Jacobi solver element-wise
   void jacobi_solver_element(T* b, T* output, int maxIter);

   // Functions for Jacobi solver matrix
   void jacobi_solver_matrix(double* b, double* output, int maxIter);
   void jacobi_decomposition(Matrix<T>* D, Matrix<T>* N);
	
   void LUDecomp(Matrix<T>& L, Matrix<T>& U);

   // Explicitly using the C++11 nullptr here
   T *values = nullptr;   
   int rows = -1;
   int cols = -1;

// We want our subclass to know about this
protected:
   bool preallocated = false;

// Private variables - there is no need for other classes 
// to know about these variables
private:

   int size_of_values = -1;   
};
