#include <iostream>
#include <math.h>
#include "Matrix.h"

// Constructor - using an initialisation list here
template <class T>
Matrix<T>::Matrix(int rows, int cols, bool preallocate): rows(rows), cols(cols), size_of_values(rows * cols), preallocated(preallocate)
{
   // If we want to handle memory ourselves
   if (this->preallocated)
   {
      // Must remember to delete this in the destructor
      this->values = new T[size_of_values];
   }
}

// Constructor - now just setting the value of our T pointer
template <class T>
Matrix<T>::Matrix(int rows, int cols, T *values_ptr): rows(rows), cols(cols), size_of_values(rows * cols), values(values_ptr)
{}

// destructor
template <class T>
Matrix<T>::~Matrix()
{
   // Delete the values array
   if (this->preallocated){
      delete[] this->values;
   }
}

// Just print out the values in our values array
template <class T>
void Matrix<T>::printValues() 
{ 
   std::cout << "Printing values" << std::endl;
	for (int i = 0; i< this->size_of_values; i++)
   {
      std::cout << this->values[i] << " ";
   }
   std::cout << std::endl;
}

// Explicitly print out the values in values array as if they are a matrix
template <class T>
void Matrix<T>::printMatrix() 
{ 
   std::cout << "Printing matrix" << std::endl;
   for (int j = 0; j< this->rows; j++)
   {  
      std::cout << std::endl;
      for (int i = 0; i< this->cols; i++)
      {
         // We have explicitly used a row-major ordering here
         std::cout << this->values[i + j * this->cols] << " ";
      }
   }
   std::cout << std::endl;
}

// Do matrix matrix multiplication
// output = this * mat_right
template <class T>
void Matrix<T>::matMatMult(Matrix& mat_left, Matrix& output)
{

   // Check our dimensions match
   if (this->cols != output.cols)
   {
      std::cerr << "Input dimensions for matrices don't match" << std::endl;
      return;
   }

   // Check if our output matrix has had space allocated to it
   if (output.values != nullptr) 
   {
      // Check our dimensions match
      if (this->rows != mat_left.cols || mat_left.rows != output.rows)
      {
         std::cerr << "Input dimensions for matrices don't match" << std::endl;
         return;
      }      
   }
   // The output hasn't been preallocated, so we are going to do that
   else
   {
      output.values = new T[this->rows * mat_left.cols];
      // Don't forget to set preallocate to true now it is protected
      output.preallocated = true;
   }

   // Set values to zero before hand
   for (int i = 0; i < output.size_of_values; i++)
   {
      output.values[i] = 0;
   }

   // Now we can do our matrix-matrix multiplication
   // CHANGE THIS FOR LOOP ORDERING AROUND
   // AND CHECK THE TIME SPENT
   // Does the ordering matter for performance. Why??
   for(int i = 0; i < this->rows; i++)
   {
      for(int k = 0; k < this->cols; k++)
      {
         for(int j = 0; j < mat_left.cols; j++)
         {            
               output.values[i * output.cols + j] += this->values[i * this->cols + k] * mat_left.values[k * mat_left.cols + j];
         }
      }
   }
}




template <class T>
void Matrix<T>::matVecMult(T* vec, T* output) {
    // Matrix vector prodcut M * b = c
    // input vec is an array which is vector b
    // output in a array of size matrix->rows to store values
     
    // Set values to zero
    for (int i = 0; i < this->rows; i++)
    {
        output[i] = 0;
    }

    // loop over rows and cols
    for (int i = 0; i < this->rows; i++) {
        for (int j = 0; j < this->cols; j++) {
            // add the right values to the output array
            output[i] += this->values[i * this->cols + j] * vec[j];
        }
    }
}


template <class T>
void Matrix<T>::vecVecsubtract(T* vec_a, T* vec_b, T* output) {
    // Vector vector subtraction Vec_A - Vec_b = output
    // all three vectors/arrays need to be same size
  
    for (int i = 0; i < this->rows; i++) {

        output[i] = vec_a[i] - vec_b[i];

    }

}


template <class T>
float Matrix<T>::RMS_norm_diff(T* vec_a, T* vec_b) {
    // RMS norm of the different of two vectors 
    // all input vectors/arrays need to be same size

    
    float sum_a = 0;

    // loop over all values in arraz
    for (int i = 0; i < this->rows; i++) {

        // add the squared difference to sum_a
        sum_a += (vec_a[i] - vec_b[i])*(vec_a[i] - vec_b[i]);
    }

    // return RMS norm of the squared difference
    return sqrt(sum_a / this->rows);




}


template <class T>
void Matrix<T>::jacobi_solver_element(T* b, T* output, int maxIter) {
    /*
    Jacobi solver using element-wise calcualtions
    Solves a linear system of equations A*x=b using an ittertive apporach
    Input:
        <T>array[]* b: RHS of the linear system
        <T>array[]* output: array in which the solution will be stored in
        int maxIter: maxiumum itterations 
    Output:
        none

    A needs to be a SPD matrix wiht no zeros on main diagonal 
    and the linear system needs to have a solution.
    Tolerance of the solver can be changed in the tol variable below

    */


    // create variables 
    double conve = 10;   // store RMS difference between x_{k} and x_{k+1}
    T* pout2 = new T[this->rows];   // store x_{k+1}
    float sum = 0; 
    int n = 0;
    //set solution tolerance to e-8
    double tol = 1.e-10;

    // initialise x_{k} 
    for (int i = 0; i < this->rows; i++)
    {
        output[i] = 1032;
    }

    // start iteration, only do maxIter steps 
    for( int n = 0; n <maxIter; n++) {
        // loop over rows
        for (int i = 0; i < this->rows; i++) {
            // set variable to zero such that it can be added to 
            sum = 0; 
            // loop over cols of matrix
            for (int j = 0; j < this->cols; j++) {
                // if i = j dont do anything because that is row value that is calcualted
                // for i not equal to j, mutiply both values and add to sum
                if (i != j) {
                    sum += this->values[i * this->cols + j] * output[j];
                }

            }
            // put the sum value in the right spot in the array with some additonal calculations
            pout2[i] = (1 / (this->values[i * this->rows + i])) * (b[i] - sum);
        }
        // swap the pointers from the x_k and x_{k+1} array such that the calculations above
        // can be repated without copying the arrays 
        std::swap(pout2, output);
        // find RMS norm of output-pout2 array 
        conve = RMS_norm_diff(pout2, output);
      
        // if rms norm is smaller than tolerance -> break loop
        if (conve < tol) {
            break;
        }
    }

    delete[] pout2;
}


template <class T>
void  Matrix<T>::jacobi_decomposition(Matrix<T>* D, Matrix<T>* N) {
    /*This function will decompose this->matirx (A) values into two matrices D and N. A = D + N
    * Matirx D will store the inverse of all values on the diagonal i.e. 1/(A[i, i])
    * Matirx N will store the lower and upper triangular parts of A
    *   Input:
        Matrix<T>* D: matrix of same dimensions as A, to store the diagonal values
        Matrix<T>* N: matrix of same dimensions as A, to store the lower and upper triangular
    Output:
        none
    */

    for (int i = 0; i < this->rows; i++) {
        for (int j = 0; j < this->cols; j++) {
            if (i == j) {
                D->values[i * this->cols + j] = (1/this->values[i * this->cols + j]);
                N->values[i * this->cols + j] = 0;
            }
            else {
                N->values[i * this->cols + j] = this->values[i * this->cols + j];
                D->values[i * this->cols + j] = 0;
            }

        }
    }
}


template <class T>
void  Matrix<T>::jacobi_solver_matrix(double* b, double* xk1, int maxIter) {
    /*
    Jacobi solver using matrix manipulation
    Solves a linear system of equations A*x=b using an ittertive apporach, where every itteration 
    mutiple matrix operations are performed: xk2 = D*(b - N * xk1)
    Input:
        <T>array[]* b: RHS of the linear system
        <T>array[]* xk1: array in which the solution will be stored in
        int maxIter: maxiumum itterations
    Output:
        none

    A needs to be a SPD matrix wiht no zeros on main diagonal
    and the linear system needs to have a solution.
    Tolerance of the solver can be changed in the tol variable below

    */

    auto* N = new Matrix<T>(this->rows, this->cols, true);
    auto* D = new Matrix<T>(this->rows, this->cols, true);



    T* xk2 = new T[this->rows];   // store x_{k+1}
    T* pvecvecarray = new T[this->rows];   // store vecVecsubtract
    T* pmatvecarray = new T[this->rows];   // store matVecMult
    //set solution tolerance to e-8
    double tol = 1.e-10;
  
    // initialize conve varible for first itteration 
    double conve = 12;


    this->jacobi_decomposition(D, N);

    
    // init start values with random values
    for (int i = 0; i < this->rows; i++)
    {
        xk1[i] = rand() % 50 + 5;
    }

    // x_{k+1} = D(b - N * x_k)

    for (int n = 0; n < maxIter; n++) {
        N->matVecMult(xk1, pmatvecarray); // K = N * x_n
        N->vecVecsubtract(b, pmatvecarray, pvecvecarray); // R = b - K
        D->matVecMult(pvecvecarray, xk2); // X_{n+1} = D * R

        // swap two pointers around
        std::swap(xk2, xk1);
        // find RMS norm of output-pout2 array 
        conve = RMS_norm_diff(xk2, xk1);

        // if rms norm is smaller than tolerance -> break loop
        if (conve < tol) {
            break;
        }
    }
  


    delete N;
    delete D;
    delete[] pmatvecarray;
    delete[] xk2;
    delete[] pvecvecarray;

}


template <class T>
void Matrix<T>::LUDecomp(Matrix& L, Matrix& U)
{

   // Check our dimensions match
   if (this->cols != L.cols || this->cols != U.cols || this->rows != L.rows || this->rows != U.rows )
   {
      std::cerr << "L and U must be of the same size as the matrix you want to decompose" << std::endl;
      return;
   }

   // Check if our L matrix has had space allocated to it
   if (L.values != nullptr) 
   {
      L.values = new T[this->rows * this->cols];
      // Don't forget to set preallocate to true now it is protected
      L.preallocated = true;
   }

    // Check if our U matrix has had space allocated to it
   if (U.values != nullptr) 
   {
      U.values = new T[this->rows * this->cols];
      // Don't forget to set preallocate to true now it is protected
      U.preallocated = true;
   }

   // Setting L = I as a starting point
   for (int i = 0; i < this->cols; i++)
   {
      for (int j = 0; j < this->rows; j++) {
         if (i == j) {
            L.values[j*this->rows + i] = 1;
         } else {
            L.values[j*this->rows + i] = 0;
         }
      }
   }

   // U = A as a starting point
   for (int i=0;i<this->size_of_values; i++)
   {
      U.values[i] = this->values[i];
   }

   for (int step = 0; step<U.rows; step++)
   {
      double factor = 0;
      for(int i = step+1; i < U.rows; i++)
      {
         factor = U.values[(i)*U.cols+step]/U.values[step*U.cols+step];
         for(int j = step; j < U.cols; j++)
         {
            U.values[i*U.cols+j] -= U.values[step*U.cols+j]*factor;
         }
         L.values[i*cols+step] = factor;
      }
   }
}
