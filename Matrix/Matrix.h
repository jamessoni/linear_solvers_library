#pragma once

template <class T>
class Matrix
{
public:

   // constructor where we want to preallocate ourselves
   Matrix(int rows, int cols, bool preallocate);
   // constructor where we already have allocated memory outside
   Matrix(int rows, int cols, T *values_ptr);
   //constructor to create SPD matrix
   //Matrix(int rows, int cols, int diag_max, int diag_min, int non_diag_max, int non_diag_min);
   Matrix(int rows, int cols, int diag_max, int diag_min);

   // destructor
   virtual ~Matrix();

   // Print out the values in our matrix
   void printValues();
   virtual void printMatrix();
   void transpose();
   // Perform some operations with our matrix
   void matMatMult(Matrix<T>& mat_right, Matrix<T>& output);
  
   bool SPDMatrixcheck();

   //Gauss-seidel solver
   void gauss_seidel(Matrix<T>& a, T* b, T* x, float tol);

   void matVecMult(T* vec, T* output);
   void vecVecsubtract(T* vec_a, T* vec_b, T* output);
   float RMS_norm_diff(T* vec_a, T* vec_b);

   // Jacobi solver element-wise
   void jacobi_solver_element(Matrix<T>* A,T* b, T* output, int maxIter, bool initialised, float tol);

   // Functions for Jacobi solver matrix
   void jacobi_solver_matrix(Matrix<T>* A, double* b, double* output, int maxIter, bool initialised, float tol);
   void jacobi_decomposition(Matrix<T>* D, Matrix<T>* N);
	
   void LUDecomp(Matrix<T>& L, Matrix<T>& U);
   void SLUDecomp(Matrix<T>* LU);
   void IPLUDecomp();
   
   void fsubstitution(Matrix<T>& L, T* y,T* b);
   void fsubstitutionLU(Matrix<T>& L, T* y,T* b);
   void bsubstitution(Matrix<T>& U, T* x, T* y);

   void LUSolve(double* b, double* output, bool inplace);
   void conjugate_gradient(T* b, T* x, int maxIter, float tol);
  
   void CholeskyDecomp(Matrix<T>* L);
   void CholeskySolve(T* b, T* x);
  
   void daxpy(int n, double alpha, double* dx, int incx, double* dy, int incy);
   void daxpytx(int n, double alpha, double* dx, int incx, double* dy, int incy);
   void dcopy(int n, double* dx, int incx, double* dy, int incy);


   // Explicitly using the C++11 nullptr here
   T *values = nullptr;   
   int rows = -1;
   int cols = -1;
   int diag_max = -1;
   int diag_min = -1;
   int non_diag_max = -1;
   int non_diag_min = -1;

// We want our subclass to know about this
protected:
   bool preallocated = false;

// Private variables - there is no need for other classes 
// to know about these variables
private:

   int size_of_values = -1;
};