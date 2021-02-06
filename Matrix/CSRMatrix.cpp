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
void CSRMatrix<T>::jacobi_solver_sparse(CSRMatrix<T>* A, T* b, T* output, int maxIter, bool initialised, float tol) {
    /*
    Jacobi solver using a CSR sparse matrix
    Solves a linear system of equations A*x=b using an iterative apporach
    Input:
        CSRMatrix<T>* A: A matrix of linear system
        <T>array[]* b: RHS of the linear system
        <T>array[]* output: array in which the solution will be stored in
        int maxIter: maxiumum itterations
        bool initialised: If true, output array was initialised with fist guess
                          If false, output array will be filled with random values
    Output:
        none
    A needs to be a SPD matrix wiht no zeros on main diagonal
    and the linear system needs to have a solution.
    */

    // create variables 
    double conve = 10;   // store RMS difference between x_{k} and x_{k+1}
    T* pout2 = new T[A->rows];   // store x_{k+1}
    float sum = 0;
    double a_ii = 0;
    double sum_RMS = 0;
    double RMS = 12;
    double* answer_check = new double[A->rows];
    //set solution tolerance to e-10


    // if not initialised fist input than use random numbers to 
    // initialise x_{k}
    if (initialised == false) {
        // initialise starting condition x_{k} 
        for (int i = 0; i < A->rows; i++)
        {
            output[i] = rand() % 500 + 50;
        }
    }

    // start iteration, only do maxIter steps 
    for (int n = 0; n < maxIter; n++) {
        // loop over rows
        for (int i = 0; i < A->rows; i++) {
            // set variable to zero such that it can be added to 
            sum = 0;
            // loop over non zero values in row
            for (int j = A->row_position[i]; j < A->row_position[i + 1]; j++) {
                // if i = j dont do anything because that is row value that is calcualted
                // for i not equal to j, mutiply both values and add to sum
                if (A->col_index[j] != i) {
                    sum += A->values[j] * output[A->col_index[j]];
                }
                if (A->col_index[j] == i) {
                    a_ii = A->values[j];
                }
                // increase j index 
            }
            // put the sum value in the right spot in the array with some additonal calculations
            pout2[i] = (1 / (a_ii)) * (b[i] - sum);
        }
        // set RMS norm summation varible to zero
        sum_RMS = 0;
        // loop over all values in arraz
        for (int i = 0; i < A->rows; i++) {

            // copy values into new array for next itteration
            output[i] = pout2[i];

        }

        A->matVecMult(output, answer_check);
        RMS = A->RMS_norm_diff(b, answer_check);
        // if rms norm is smaller than tolerance -> break lool
        if (RMS < tol) {
            break;
        }

    }
    delete[] pout2;
    delete[] answer_check;
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
void CSRMatrix<T>::gauss_seidel_sparse(CSRMatrix<T>& a, T* b, T* x_init, float tol)
{
    /*
    Gauss-seidel solver implementation
        Solves a linear system of equations A * x = b using an iterative apporach
        Input :
            <T>array[] * a : CSR Matrix input
            <T>array[] * b : RHS of the linear system
            <T>array[] * x_init : array in which the solution will be stored in
            double tol : tolerance of the solver

        Output :
            none

        A needs to be a diagonal dominant SPD matrix with no zeros on main diagonal
        and the linear system needs to have a solution.
        For gauss-seidel to be ran on maximum iterations - uncomment '//while (iter < iter_max)' below.
        int iter_max : maximum number of iterations can be altered below

    */

    int iter_max = 500;
    int iter = 0;
    double conve = 10;
    double A_ii = 0;

    T* pout2 = new T[this->rows];

    //x_init has initialised 0's as an initial 'guess'
    for (int i = 0; i < this->rows; i++)
    {
        x_init[i] = 0;
    }

    while (conve > tol)
        //while (iter < iter_max)
    {
        iter += 1;
        //updating pout2 as previous iteration x.values
        //enables convergence parameter to be checked against predefined tolerance
        for (int i = 0; i < this->rows; i++)
        {
            pout2[i] = x_init[i];
        }
        //looping over the rows of A
        for (int i = 0; i < this->rows; i++)
        {

            //A_ii = 0;
            //looping through storing the diagonal entry in A_ii
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
            x_init[i] = b[i] / A_ii;
            //using compressed sparse row matVecMult 
            //looping over the columns 
            for (int val_index = a.row_position[i]; val_index < a.row_position[i + 1]; val_index++)
            {
                if (a.col_index[val_index] == i)
                {
                    continue;
                }
                //computing x[i] = x[i] - ((A[i][j] * x[j]) / A[i][i])
                x_init[i] = x_init[i] - (a.values[val_index] * x_init[col_index[val_index]]) / A_ii;
            }
        }
        // calculating the rms norm difference between the previous
        // and current iterations of x
        conve = RMS_norm_diff(pout2, x_init);
        //std::cout << "\nconvergence values: " << conve;
    }

    delete[] pout2;
    
  }

template <class T>
int CSRMatrix<T>::getv(int row,int col)
{
    for (int v = this->row_position[row]; v < this->row_position[row + 1]; v++) //check all nnzvs in a row
    {
        if (this->col_index[v] == col) // check if the nnzv in the row has the column index we are looking for
        {
            return v;
        }
    }
return INT_MIN;
}

// Do matrix matrix multiplication
// output = this * mat_right
template <class T>
CSRMatrix<T>* CSRMatrix<T>::matMatMult(CSRMatrix<T>& mat_right)
{
    //better dimcheck

   // Check our dimensions match
   if (this->cols != mat_right.cols)
   {
      std::cerr << "Input dimensions for matrices don't match" << std::endl;
      return;
   }

    int nnzs=0;
    vector<T> values;
    vector<T> row_position;
    vector<T> col_index;

    row_position.push_back(0); //first 0
    T product = 0;
    int v2;
	    
        //looping through each row of the left matrix
       for (int i = 0; i < this->rows; i++)
       {
		   //loop through each column of the right matrix
           for (int j = 0; j < mat_right.cols; j++)
           {
               product = 0;
			   //if there are two values in row 3, we go 4 to 6
                for (int v = this->row_position[i]; v < this->row_position[i + 1]; v++)
                {
					// //for each nnz found check if we have an nnz in the right place in the right matrix too
                    v2 = mat_right.getv(this->col_index[v],  j);
                    if (v2 > 0)
                    {
                        product = product + this->values[v] * mat_right.values[v2];
                    }
                }
                if (product!=0)
                {
                    nnzs+=1;
                    values.push_back(product);
                    col_index.push_back(j);
                }
           }
        row_position.push_back(nnzs);
        }
    auto* toreturn = new CSRMatrix(this->rows,this->cols,nnzs,true);
    
    //really really ugly setting values and col indexes
    for (int i = 0; i<values.size();i++)
    {
        toreturn->values[i] = values[i];
        toreturn->col_index[i] = col_index[i];
    }

    //really ugly setting row position values
    for (int i = 0; i<row_position.size();i++)
    {
        toreturn->row_position[i] = row_position[i];
    }

    return toreturn;
}

template <class T>
void CSRMatrix<T>::CholeskyDecomp()
{
    int nnzs=0;
    vector<T> values;
    vector<T> row_position;
    vector<T> col_index;
    row_position.push_back(0);

//     for (int i = 0; i < this->rows; i++)
//    {  
//       //all the non diagonal elements first
//         for (int v = this->row_position[i]; v < this->row_position[i + 1]; v++)
//         {
//          double sigma = 0;
//          // non-diagonals L(i, j) 
//          for (int p = 0; p < j; p++)
//          {
//             sigma += L->values[i*this->cols + p] * L->values[j*this->cols + p]; 
//          }
//          L->values[i*this->cols + j] = ((this->values[i*this->cols + j] - sigma) / L->values[j*this->cols + j]);
//       }

//       //Then the diagonal element
//       double sigma = 0;
//       for (int p = 0; p < i; p++) 
//       {
//          sigma += L->values[i*this->cols + p] * L->values[i*this->cols + p];
//       }
      
//       L->values[i*this->cols + i] = sqrt(this->values[i*this->cols + i] - sigma);

//       for (int j = i+1; j < this->cols; j++)
//       {
//          L->values[i*this->cols + j] = 0;
//       }
//    }
}