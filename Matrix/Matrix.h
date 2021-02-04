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
    
   
   void jacobi_solver_element(double* b, float* output, int maxIter, float tol);

   void jacobi_solver_matrix(double* b, double* output, int maxIter, float tol);
   void jacobi_decomposition(Matrix<T>* D, Matrix<T>* N);

   void LUDecomp(Matrix<T>& L, Matrix<T>& U);
   void SLUDecomp(Matrix<T>& LU);
   void IPLUDecomp();
   
   void fsubstitution(Matrix<T>& L, T* y,T* b);
   void bsubstitution(Matrix<T>& U, T* x, T* y);

   void LUSolve(double* b, double* output, bool inplace);
   void conjugate_gradient(T* b, T* x, int maxIter, float tol);

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
