#include <iostream>
#include "Matrix.h"

using namespace std;

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

// Matrix vector prodcut M * b = c
// input vec is an array which is vector b
// output in a array of size matrix->rows to store values
template <class T>
void Matrix<T>::matVecMult(T* vec, T* output) {
     
    // Set values to zero
    for (int i = 0; i < this->rows; i++)
    {
        output[i] = 0;
    }

    for (int i = 0; i < this->rows; i++) {
        for (int j = 0; j < this->cols; j++) {
            output[i] += this->values[i * this->cols + j] * vec[j];
        }
    }
}

template <class T>
void Matrix<T>::vecVecsubtract(T* vec_a, T* vec_b, T* output) {

    if ((sizeof(vec_a) / sizeof(*vec_a)) != (sizeof(vec_b) / sizeof(*vec_b))) {
        cout << "vector are not the same size for subtraction";
        return;
    }


    for (int i = 0; i < this->rows; i++) {

        output[i] = vec_a[i] - vec_b[i];

    }

}


template <class T>
void Matrix<T>::jacobi_solver_element(double* b, float* output, int maxIter, float tol) {

    const int rows_c = (sizeof(b) / sizeof(*b));
    float conve = 10;
    float output2[rows_c];
    float* pout2 = output2;
    float sum = 0;
    


    for (int i = 0; i < this->rows; i++)
    {
        output[i] = 1;
    }

    cout << "\nwewew";

    

    while (conve != tol) {
        for (int n = 0; n < maxIter; n++) {
            for (int i = 0; i < this->rows; i++) {
                sum = 0;
                for (int j = 0; j < this->cols; j++) {
                        if (i != j) {
                            sum += this->values[i * this->cols + j] * output[j];
                        }

                }

                pout2[i] = (1 / (this->values[i * this->cols + i])) * (sum + b[i]);
            }


            cout << "output  " << pout2 << "\n";
            cout << "output2  " << output2 << "\n";
            
                      
            swap(pout2, output);
            
           


        }




    }
 
}


template <class T>
void  Matrix<T>::jacobi_decomposition(Matrix<T>* D, Matrix<T>* N) {

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
void  Matrix<T>::jacobi_solver_matrix(double* b, double* xk1, int maxIter, float tol) {

    auto* N = new Matrix<T>(this->rows, this->cols, true);
    auto* D = new Matrix<T>(this->rows, this->cols, true);
    const int rows_c = (sizeof(b) / sizeof(*b));
    T output2[rows_c];
    T* xk2 = output2;
    T vecvecarray[rows_c];
    T* pvecvecarray = vecvecarray;
    T matvecarray[rows_c];
    T* pmatvecarray = matvecarray;
    T conve = 12;

    cout << "matrix jacobi";

    this->jacobi_decomposition(D, N);

    
    // init start values
    for (int i = 0; i < this->rows; i++)
    {
        xk1[i] = 1;
    }

    // x_{k+1} = D(b - N * x_k)
    while (conve != tol) {
        for (int n = 0; n < maxIter; n++) {
            N->matVecMult(xk1, pmatvecarray); // O = N * x_n
            N->vecVecsubtract(b, pmatvecarray, pvecvecarray); // R = b - O
            N->matVecMult(pvecvecarray, xk2); // X_{n+1} = D * R
        }
    }


    delete N;
    delete D;

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

template <class T>
void Matrix<T>::LUSolve(double* b, double* output)
{
   auto *L = new Matrix<double>(this->rows, this->cols, true);
   auto *U = new Matrix<double>(this->rows, this->cols, true);
   this->LUDecomp(*L,*U);

   //std::cout<<b[0]<<b[1]<<b[2];
   int size = 3;
   double y[3];
   for (int i = 0; i<size;i++)
   {
      double sum =0;
      for (int j = 0; j<size; j++)
      {
         sum+= L->values[i*L->cols + j]*y[j];
      }
      y[i]=(b[i]-sum)/L->values[i*(L->cols+1)];
   }

   for (int i = size-1; i>=0;i--)
   {
      double sum =0;
      for (int j = 0; j<size; j++)
      {
         sum+= U->values[i*U->cols + j]*output[j];
      }

      output[i]=(y[i]-sum)/U->values[i*(U->cols+1)];
   }

   delete L;
   delete U;
}