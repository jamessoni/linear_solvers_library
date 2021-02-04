#include <iostream>
#include "CSRMatrix.h"

// Constructor - using an initialisation list here
template <class T>
CSRMatrix<T>::CSRMatrix(int rows, int cols, int nnzs, bool preallocate): Matrix<T>(rows, cols, false), nnzs(nnzs)
{
   // If we don't pass false in the initialisation list base constructor, it would allocate values to be of size
   // rows * cols in our base matrix class
   // So then we need to set it to the real value we had passed in
   this->preallocated = preallocate;

   // If we want to handle memory ourselves
   if (this->preallocated)
   {
      // Must remember to delete this in the destructor
      this->values = new T[this->nnzs];
      this->row_position = new int[this->rows+1];
      this->col_index = new int[this->nnzs];
   }
}

// Constructor - now just setting the value of our T pointer
template <class T>
CSRMatrix<T>::CSRMatrix(int rows, int cols, int nnzs, T *values_ptr, int *row_position, int *col_index): Matrix<T>(rows, cols, values_ptr), nnzs(nnzs), row_position(row_position), col_index(col_index)
{}

// destructor
template <class T>
CSRMatrix<T>::~CSRMatrix()
{
   // Delete the values array
   if (this->preallocated){
      delete[] this->row_position;
      delete[] this->col_index;
   }
   // The super destructor is called after we finish here
   // This will delete this->values if preallocated is true
}

// Explicitly print out the values in values array as if they are a matrix
template <class T>
void CSRMatrix<T>::printMatrix() 
{ 
   std::cout << "Printing matrix" << std::endl;
   std::cout << "Values: ";
   for (int j = 0; j< this->nnzs; j++)
   {  
      std::cout << this->values[j] << " ";      
   }
   std::cout << std::endl;
   std::cout << "row_position: ";
   for (int j = 0; j< this->rows+1; j++)
   {  
      std::cout << this->row_position[j] << " ";      
   }
   std::cout << std::endl;   
   std::cout << "col_index: ";
   for (int j = 0; j< this->nnzs; j++)
   {  
      std::cout << this->col_index[j] << " ";      
   }
   std::cout << std::endl;   
}

// Do a matrix-vector product
// output = this * input
template<class T>
void CSRMatrix<T>::matVecMult(T *input, T *output)
{
   if (input == nullptr || output == nullptr)
   {
      std::cerr << "Input or output haven't been created" << std::endl;
      return;
   }

   // Set the output to zero
   for (int i = 0; i < this->rows; i++)
   {
      output[i] = 0.0;
   }

   // Loop over each row
   for (int i = 0; i < this->rows; i++)
   {
      // Loop over all the entries in this col
      for (int val_index = this->row_position[i]; val_index < this->row_position[i+1]; val_index++)
      {
         // This is an example of indirect addressing
         // Can make it harder for the compiler to vectorise!
         output[i] += this->values[val_index] * input[this->col_index[val_index]];

      }
   }
}


// Do matrix matrix multiplication
// output = mat_left * this
template <class T>
void CSRMatrix<T>::matMatMult(CSRMatrix<T>& mat_left, CSRMatrix<T>& output)
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
      std::cerr << "OUTPUT HASN'T BEEN ALLOCATED" << std::endl;

   }

   // HOW DO WE SET THE SPARSITY OF OUR OUTPUT MATRIX HERE??
}



template<class T>
void CSRMatrix<T>::dense2sparse(Matrix<T>& tosparsify, CSRMatrix<T>* output)
{
    int nnzvs = 0;

    for (int i = 0; i < tosparsify.rows * tosparsify.cols; i++)
    {
        if (!tosparsify.values[i] == 0)
        {
            nnzvs++;
        }
    }


    int nnzv = 0;


    T* vals = new T[nnzvs];
    int* col_index = new int[nnzvs];
    int* row_position = new int[nnzvs + 1];
    row_position[0] = 0;

    for (int i = 0; i < tosparsify.rows; i++)
    {
        for (int j = 0; j < tosparsify.cols; j++)
        {
            if (!tosparsify.values[i * tosparsify.cols + j] == 0)
            {
                vals[nnzv] = tosparsify.values[i * tosparsify.cols + j];
                col_index[nnzv] = (i * tosparsify.cols + j) % tosparsify.cols;
                nnzv++;
            }
        }
        row_position[i + 1] = nnzv;
    }


    output->nnzs = nnzvs;
    output->values = vals;
    output->col_index = col_index;
    output->row_position = row_position;

    // delete to pointers to memory location but not the array in the memory

}


template <class T>
void CSRMatrix<T>::jacobi_solver_sparse(T* b, T* output, int maxIter, bool initialised) {
    /*
    Jacobi solver using element-wise calcualtions
    Solves a linear system of equations A*x=b using an ittertive apporach
    Input:
        <T>array[]* b: RHS of the linear system
        <T>array[]* output: array in which the solution will be stored in
        int maxIter: maxiumum itterations
        bool initialised: If true, output array was initialised with fist guess
                          If false, output array will be filled with random values
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
    int j_index = 0;
    double a_ii = 0;
    double sum_RMS = 0;
    //set solution tolerance to e-10
    double tol = 1.e-10;


    // if not initialised fist input than use random numbers to 
    // initialise x_{k}
    if (initialised == false) {
        // initialise starting condition x_{k} 
        for (int i = 0; i < this->rows; i++)
        {
            output[i] = rand() % 500 + 50;
        }
    }


    // start iteration, only do maxIter steps 
    for (int n = 0; n < maxIter; n++) {
        // loop over rows
        for (int i = 0; i < this->rows; i++) {
            // set variable to zero such that it can be added to 
            sum = 0;
            j_index = 0;
            // loop over non zero values in row
            for (int j = this->row_position[i]; j < this->row_position[i + 1]; j++) {
                // if i = j dont do anything because that is row value that is calcualted
                // for i not equal to j, mutiply both values and add to sum
                if (this->col_index[j] != i) {
                    //cout << "this->values[j] " << this->values[j] << "  j  " << j_index << " output[j] " << output[j_index] << endl;
                    sum += this->values[j] * output[j_index];
                }

                if (this->col_index[j] == i) {
                    a_ii = this->values[j];
                }

                j_index++;
            }
            // put the sum value in the right spot in the array with some additonal calculations
            //cout << "b[i] " << b[i] << "  sum " << sum << endl;
            pout2[i] = (1 / (a_ii)) * (b[i] - sum);
            //cout << " pout2 cal " << pout2[i] << endl;
        }

        // set RMS norm summation varible to zero
        sum_RMS = 0;

        // loop over all values in arraz
        for (int i = 0; i < this->rows; i++) {

            //cout << pout2[i] << "  " << output[i] << endl;
            // add the squared difference to sum_a for the RMS convergence between the
            // different itterations
            sum_RMS += (pout2[i] - output[i]) * (pout2[i] - output[i]);

            // copy values into new array for next itteration
            output[i] = pout2[i];


        }

        // RMS norm of the squared difference
        conve = sqrt(sum_RMS / this->rows);

        //cout << "N " << n << " conve : " << conve << endl;

        // if rms norm is smaller than tolerance -> break loop
        if (conve < tol) {
            break;
        }
    }

    delete[] pout2;
}

template <class T>
float CSRMatrix<T>::RMS_norm_diff(T* vec_a, T* vec_b)
{
    // RMS norm of the different of two vectors 
    // all input vectors/arrays need to be same size

    float sum_a = 0;

    // loop over all values in arraz
    for (int i = 0; i < this->rows; i++)
    {
        // add the squared difference to sum_a
        sum_a += (vec_a[i] - vec_b[i]) * (vec_a[i] - vec_b[i]);

    }

    // return RMS norm of the squared difference
    return sqrt(sum_a / this->rows);
}

template <class T>
void CSRMatrix<T>::gauss_seidel(CSRMatrix<T>& a, Matrix<T>& b, Matrix<T>& x_init)
{
    //   Gauss-seidel implementation
    //   Method for solving a linear system, Ax = b, where A is a positive definite matrix (sparse matrix)
    //   Both convergence tolerance and fixed iteration methodologies presented 
    //   as convergence criteria
    //   Convergence tolerance and iteration number are predefined in the variables tol and iter_max below.

    double tol = 1e-10;
    int iter_max = 500;
    int iter = 0;
    double conve = 10;
    double A_ii = 0;

    T* pout2 = new T[x_init.rows];

    //x_init has initialised 0's
    for (int i = 0; i < x_init.rows; i++)
    {
        x_init.values[i] = 0;
    }

    while (conve > tol)
        //while (iter < iter_max)
    {
        iter += 1;
        std::cout << "\niteration: " << iter;
        //updating pout2 as previous iteration x.values
        //enables convergence parameter to be checked against predefined tolerance
        for (int i = 0; i < x_init.rows; i++)
        {
            pout2[i] = x_init.values[i];
        }

        for (int i = 0; i < this->rows; i++)
        {

            //A_ii = 0;
            //looping through diagonal entries for A_ii
            for (int val_index = a.row_position[i]; val_index < a.row_position[i + 1]; val_index++)
            {

                if (a.col_index[val_index] == i)
                {
                    A_ii = a.values[val_index];
                    //once A_ii its value stored - break to prevent overwriting
                    if (A_ii != 0)
                    {
                        break;
                    }
                }
            }
            x_init.values[i] = b.values[i] / A_ii;
            //using compressed sparse row matVecMult 
            for (int val_index = a.row_position[i]; val_index < a.row_position[i + 1]; val_index++)
            {
                if (a.col_index[val_index] == i)
                {
                    continue;
                }
                x_init.values[i] = x_init.values[i] - (a.values[val_index] * x_init.values[col_index[val_index]]) / A_ii;
            }
        }
        conve = RMS_norm_diff(pout2, x_init.values);
        std::cout << "\nconvergence values: " << conve;
    }
    std::cout << std::endl;

    delete[] pout2;

}