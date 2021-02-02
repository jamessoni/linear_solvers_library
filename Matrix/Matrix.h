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

   // Perform some operations with our matrix
   void matMatMult(Matrix<T>& mat_left, Matrix<T>& output);
   double RMS_norm_diff(T* vec_a, T* vec_b);
   bool SPDMatrixcheck();

   //Gauss-seidel solver
   void gauss_seidel(Matrix<T>& a, Matrix<T>& b, Matrix<T>& x);

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