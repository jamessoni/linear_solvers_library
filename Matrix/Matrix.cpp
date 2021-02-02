#include <iostream>
#include "Matrix.h"
#include <vector>

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

////Constructor - creating the SPD matrix
template <class T>
Matrix<T>::Matrix(int rows, int cols, int diag_max, int diag_min) : rows(rows), cols(cols),
diag_max(diag_max), diag_min(diag_min), non_diag_max(non_diag_max), non_diag_min(non_diag_min), size_of_values(rows* cols), preallocated(true)
{
    this->values = new T[this->rows * this->cols];
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            if (i == j)
            {
                this->values[i * cols + j] = (rand() % diag_max) + diag_min;
            }
            else
            {
                this->values[i * cols + j] = 1;
            }
        }
    }
}

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
bool Matrix<T>::SPDMatrixcheck()
{
    //this method returns true if the matrix is symmetric 
    // AND the sum of the non diagonal entries of the row
    // is less than the smallest diagonal entry

    int sum_nondiag{};
    int diag_val{};
    bool SPD{ false };
    //bool weak_test{ false };

    //testing for symmetry in matrix
    for (int i = 0; i < rows; i++){
        for (int j = 0; j < cols; j++){
            if (i == j){
                continue;
            }
            else{
                if (this->values[i * rows + j] == this->values[j * cols + i]){
                    SPD = true;
                }
                else{
                    SPD = false;
                }
            }
        }
    }
    //SPD test (weak) that enables convergence to occur
    for (int i = 0; i < rows; i++){
        for (int j = 0; j < cols; j++){
            if (i == j){
                diag_val = this->values[i * rows + j];
                continue;
            }
            else{
                sum_nondiag += this->values[i * rows + j];
            }
        }
    }
    if (diag_val <= sum_nondiag) { 
        SPD = false; 
    }
    if (SPD == true){
        std::cout << "\nPasses weak SPD check";
}
    else {
        std::cout << "\nMay not converge with chosen solvers";
    }
    return SPD;
}

template <class T>
double Matrix<T>::RMS_norm_diff(T* vec_a, T* vec_b) 
{
    // RMS norm of the different of two vectors 
    // all input vectors/arrays need to be same size

    float sum_a = 0;

    // loop over all values in arraz
    for (int i = 0; i < this->rows; i++) {

        // add the squared difference to sum_a
        sum_a += (vec_a[i] - vec_b[i]) * (vec_a[i] - vec_b[i]);
    }

    // return RMS norm of the squared difference
    return sqrt(sum_a / this->rows);

}

template <class T>
void Matrix<T>::gauss_seidel(Matrix<T>& a, Matrix<T>& b, Matrix<T>& x_init)
{
    //Both convergence tolerance and fixed iteration methodologies presented
    //design decision that although set iteration gives a good method for
    // large matrices. For small matrices it may not be the most effective

    double tol = 1e-5;
    int iter_max = 500;
    int iter = 0;
    double conve = 10;

    T* pout2 = new T[x_init.rows];

    std::shared_ptr<Matrix<T>> y(new Matrix<T>(x_init.rows, x_init.cols, true));
    //filling y
     for (int i = 0; i < y->rows * y->cols; i++)
    {
        y->values[i] = 0;
    }

    // std::shared_ptr<Matrix<T>> x(new Matrix<T>(x_init.rows, x_init.cols, true));

     while (conve > tol)
     {
         iter += 1;
         std::cout << "\niteration: " << iter;
         for (int i = 0; i < x_init.rows; i++)
         {
             pout2[i] = x_init.values[i];
         }

         for (int i = 0; i < this->rows; i++)
         {
             y->values[i] = b.values[i] / a.values[i * rows + i];
             for (int j = 0; j < this->cols; j++)
             {
                 if (i == j)
                 {
                     continue;
                 }
                 y->values[i] = y->values[i] - ((a.values[i * rows + j] / a.values[i * rows + i]) * x_init.values[j]);
                 x_init.values[i] = y->values[i];
             }
         }
         conve = RMS_norm_diff(pout2, x_init.values);
         std::cout << "\nconvergence values: " << conve;
     }
     std::cout << std::endl;
    
     delete[] pout2;

    // Same implementation as above which is only based off number of iterations

    //while (iter_max > 0)
    //{
    //    //for (int i = 0;x->values[i] = x_init[i];
    //    for (int i = 0; i < this->rows; i++)
    //    {
    //        y->values[i] = b.values[i] / a.values[i * rows + i];
    //        for (int j = 0; j < this->cols; j++)
    //        {
    //            if (i == j)
    //            {
    //                continue;
    //            }
    //            y->values[i] = y->values[i] - ((a.values[i * rows + j] / a.values[i * rows + i]) * x_init.values[j]);
    //            x_init.values[i] = y->values[i];
    //        }
    //    }   
    //    iter -= 1;
    //}  
    //std::cout << std::endl;

}