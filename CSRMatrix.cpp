#include <iostream>
#include "CSRMatrix.h"
#include <vector>

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
void CSRMatrix<T>::jacobi_solver_sparse(CSRMatrix<T>* A, T* b, T* x, int maxIter, bool initialised, float tol) {
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
    // Check our dimensions match
    if (this->cols != this->rows)
    {
        std::cerr << "Ensure dimensions (rows and columns) of Matrix A match" << std::endl;
    }
    else {
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
                x[i] = rand() % 500 + 50;
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
                        sum += A->values[j] * x[A->col_index[j]];
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
                x[i] = pout2[i];

            }

            A->matVecMult(x, answer_check);
            RMS = A->RMS_norm_diff(b, answer_check);
            // if rms norm is smaller than tolerance -> break lool
            if (RMS < tol) {
                break;
            }

        }
        delete[] pout2;
        delete[] answer_check;
    }
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
void CSRMatrix<T>::gauss_seidel_sparse(CSRMatrix<T>* A, T* b, T* x, int maxIter, float tol)
{
    /*
    Gauss-seidel solver implementation
        Solves a linear system of equations A * x = b using an iterative apporach
        Input :
            <T>array[] * a : CSR Matrix input
            <T>array[] * b : RHS of the linear system
            <T>array[] * x : array in which the solution will be stored in
            double tol : tolerance of the solver

        Output :
            none

        A needs to be a diagonal dominant SPD matrix with no zeros on main diagonal
        and the linear system needs to have a solution.
        For gauss-seidel to be ran on maximum iterations - uncomment '//while (iter < iter_max)' below.
        int iter_max : maximum number of iterations can be altered below

    */

    // Check our dimensions match
    if (this->cols != this->rows)
    {
        std::cerr << "Ensure dimensions (rows and columns) of Matrix A match" << std::endl;
    }
    else {

        int iter = 0;
        double conve = 10;
        double A_ii = 0;
        double* answer_check = new double[A->rows];
        float RMS;
        T* pout2 = new T[this->rows];

        //x_init has initialised 0's as an initial 'guess'
        for (int i = 0; i < this->rows; i++)
        {
            x[i] = 0;
        }

        int count = 0;
        while (count < maxIter)
            //while (iter < iter_max)
        {
            count += 1;
            //updating pout2 as previous iteration x.values
            //enables convergence parameter to be checked against predefined tolerance
            for (int i = 0; i < this->rows; i++)
            {
                pout2[i] = x[i];
            }
            //looping over the rows of A
            for (int i = 0; i < this->rows; i++)
            {

                //A_ii = 0;
                //looping through storing the diagonal entry in A_ii
                for (int val_index = A->row_position[i]; val_index < A->row_position[i + 1]; val_index++)
                {

                    if (A->col_index[val_index] == i)
                    {
                        A_ii = A->values[val_index];
                        //once A_ii its value stored - break to prevent overwriting
                        if (A_ii != 0)
                        {
                            break;
                        }
                    }
                }
                x[i] = b[i] / A_ii;
                //using compressed sparse row matVecMult 
                //looping over the columns 
                for (int val_index = A->row_position[i]; val_index < A->row_position[i + 1]; val_index++)
                {
                    if (A->col_index[val_index] == i)
                    {
                        continue;
                    }
                    //computing x[i] = x[i] - ((A[i][j] * x[j]) / A[i][i])
                    x[i] = x[i] - (A->values[val_index] * x[col_index[val_index]]) / A_ii;
                }
            }
            // calculating the rms norm difference between the previous
            // and current iterations of x
            A->matVecMult(x, answer_check);
            RMS = A->RMS_norm_diff(b, answer_check);
            // if rms norm is smaller than tolerance -> break loop
            if (RMS < tol) {
                break;
            }

        }

        delete[] pout2;
        delete[] answer_check;
    }
  }

 //Do matrix matrix multiplication
 //output = this * mat_right
//template <class T>
//CSRMatrix<T>* CSRMatrix<T>::matMatMult(CSRMatrix<T>& mat_right)
//{
//    //add better dimcheck
//
//   // Check if our dimensions match
//   if (this->cols != mat_right.rows)
//   {
//      std::cerr << "Input dimensions for matrices don't match" << std::endl;
//   }
//
//    //allocating some space for our CSRMatrix.
//    //row_position can be allocated as a vector because we know how many rows there are already
//    //values and col_index cannot be allocated as vectors because we do not know how many nnzs there are
//
//    int nnzs=0;
//    vector<T> values;
//    vector<int> col_index;
//    auto* row_position = new int[this->rows+1];
//    
//    //first element of row_position is always 0
//    row_position[0]=0;
//
//    //Keeping track of our row/column product
//    T product;
//    int v2;
//	    
//        //looping through every row of the left matrix
//       for (int i = 0; i < this->rows; i++)
//       {
//		   //looping through every column of the right matrix
//           for (int j = 0; j < mat_right.cols; j++)
//           {
//               product = 0;
//               //looping through all values in the left matrix's row
//                for (int v = this->row_position[i]; v < this->row_position[i + 1]; v++)
//                {
//                    int index = -1;
//					// for each nnz found check if we have an nnz in the right place in the right matrix too
//                    //check all nnzvs in the row
//                    for (int v = this->row_position[this->col_index[v]]; v < this->row_position[this->col_index[v] + 1]; v++)
//                    {
//                        // If any of them has the col index we are looking for
//                        if (this->col_index[v] == j) 
//                        {
//                            //return the value index
//                            index = v;
//                        }
//                    }
//                    //if we find it
//                    if (v2 > 0)
//                    {
//                        //multiply and add to the product
//                        product += this->values[v] * mat_right.values[v2];
//                    }
//                }
//
//                //if the total product is not 0, then add the value to our work in progress CSRMatrix
//                if (product!=0)
//                {
//                    nnzs+=1;
//                    values.push_back(product);
//                    col_index.push_back(j);
//                }
//            }
//
//            //at the end of every row loop, update row_position;
//            row_position[i+1] = nnzs;
//        }
//
//    auto* toreturn = new CSRMatrix(this->rows,this->cols,nnzs,true);
//    
//    //copying our values over to our newly created CSRmatrix
//
//    copy(values.begin(),values.end(),toreturn->values);
//    copy(col_index.begin(),col_index.end(),toreturn->col_index);
//    toreturn->row_position = row_position;
//
//    //and finally returning the CSRmatrix itself
//    return toreturn;
//}

template <class T>
void CSRMatrix<T>::transposeiflower()
{
    //CAREFUL! THIS ONLY WORKS FOR LOWER TRIANGULAR MATRICES
    //Very slow, O(n^2), sure I could do better with more time

    int nnzs = this->nnzs;
    auto* values =  new T[nnzs];
    auto* col_index = new int[nnzs];
    auto* row_position = new int[this->cols+1];

    for (int i = 0; i<this->cols+1;i++)
    {
        row_position[i] = 0;
    }

    int nnzsprocessed = 0;
    int search = 0;

    while (nnzsprocessed<nnzs)
    {
        for (int i=0; i<nnzs;i++)
        {
            if (this->col_index[i] == search)
            {
                values[nnzsprocessed] = this->values[i];
                for (int j = this->col_index[i]+1;j<this->cols+1;j++)
                {
                    row_position[j]++;
                }
                for (int row=this->rows;row>=0;row--)
                {
                    col_index[nnzsprocessed] = row;
                    if (i+1>this->row_position[row])
                    {
                        break;
                    }
                }
                nnzsprocessed++;
            }
        }
        search++;
    }
    this->col_index = col_index;
    this->values = values;
    this->row_position = row_position;
}


template <class T>
CSRMatrix<T>* CSRMatrix<T>::CholeskyDecomp()
{   
    //lots of memory taken up but I could not find any other solution
    //I briefly tried vectors, but it was just weird
    int nnzs=0;
    auto* values = new T[this->rows*this->cols];
    auto* row_position = new int[this->rows+1];
    auto* col_index = new int[this->rows * this->cols];

    //storing ALL diagonal values makes the calculation a LOT easier
    auto* diagonal = new T[this->rows];

    row_position[0] = 0; //first element is always 0
    row_position[1] = 0; //second element will also be 0 for the loop

    int this_index;

    //Looping through all rows
    for (int i  = 0; i < this->rows; i++)
    {
        //Keep track of the val index of the starting matrix
        this_index = this->row_position[i];

        //all the non diagonal elements first
        for (int j = 0; j < i; j++)
        {   
            //Lij = aij - sigma((Lip*Ljp), for p<j) // Ljj

            //Numerator
            T num = 0;

            //Find the sum of L[i][p]*L[j][p] for all p<j
            for (int v = row_position[i]; (v < row_position[i+1] && col_index[v] < j); v++) //first loop
            {
                for (int m = row_position[j]; (m < row_position[j+1] && col_index[m] <= col_index[v]); m++) //second loop
                {
                    if (col_index[m] == col_index[v]) //if there is a correspondence
                    {
                        num -= values[v] * values[m]; //negative because we need to subtract sigma from the numerator
                    }
                }
            }

            //Check if aij exists
            if (this->col_index[this_index] == j)
            {   
                //if it does add it to the numerator
                num +=this->values[this_index];
                this_index++;
            }

            //If our numerator != 0, we have a non zero value to add to our matrix
            if (num != 0)
            {
                //We divide by Ljj
                T val = num / diagonal[j];
                
                //And then we add the result
                values[nnzs] = val;
                col_index[nnzs] = j;
                row_position[i+1]++;
                nnzs++;
            }
        }

        //out of the loop for non-diagonal elements. We proceed to the diagonal element in our row
        
        //Lii = SQRT(aii - sigma((Lip^2), for p<i))
        T tosquareroot = 0;

        //looping
        for (int v = row_position[i]; v < row_position[i+1]; v++)
        {
            //squaring and detracting from tosquareroot
            tosquareroot -= pow(values[v],2);
        }

        //adding aii
        tosquareroot += this->values[this_index];

        //squarerooting to find our value
        T val = sqrt(tosquareroot);

        //Add elements to matrix
        values[nnzs] = val;
        col_index[nnzs] = i;
        row_position[i+1]++;
        row_position[i+2]=row_position[i+1];

        diagonal[i] = val;

        nnzs++;
        
    }

    //Once we have "built" our vectors, we build our matrix and return it
    auto* toreturn = new CSRMatrix(this->rows,this->cols,nnzs,true);
    
    toreturn->values = values;
    toreturn->col_index = col_index;
    toreturn->row_position = row_position;

    return toreturn;
}

template <class T>
void CSRMatrix<T>::fsubstitution(T* b, T* y)
{
    T sigma; //sum of product of before-diag elements with construction-in-progress y vector
    T diag; //Keeping track of the diagonal for efficiency

    //Looping through all rows
    for (int i = 0; i < this->rows; i++)
    {
        //setting the sum to 0
        sigma = 0;

        //looping through values in row
        for (int v = this->row_position[i]; v < this->row_position[i+1]; v++)
        {
            //column index j of the value we are currently considering
            int j = this->col_index[v]; 

            if (i==j) //if i=j, then it's a diagonal value. We store it for later use
            {
                diag = this->values[v];
            }
            if (j < i) //Stopping at the diagonal
            {
                //multiplying value in column j of the Cholesky matrix by the value of the row v in y, adding to total
                sigma += this->values[v]  * y[j];
            }
            if (j > i) //don't want to keep looping
            {
                break;
            }
        }
        y[i] = (b[i] - sigma) / diag; //updating y
    }
}

template <class T>
void CSRMatrix<T>::bsubstitution(T* y, T* x)
{
    T sigma; //sum of product of before-diag elements with construction-in-progress y vector
    T diag; //Keeping track of the diagonal for efficiency

    //Looping through all rows, bottom to top
    for (int i = this->rows - 1; i >= 0; i--)
    {
        //setting the sum to 0
        sigma = 0;

        //looping through values in row
        for (int v = this->row_position[i]; v < this->row_position[i + 1]; v++)
        {
            //column index j of the value we are currently considering
            int j = this->col_index[v];

            if (j > i)  //Stopping at the diagonal
            {
                //multiplying value in column j of the Cholesky matrix by the value of the row v in y, adding to total
                sigma += this->values[v] * x[j];
            }
            if (j == i) //if i=j, then it's a diagonal value. We store it for later use
            {
                diag = this->values[v];
            }
            if (j < i) //don't want to keep looping
            {
                break;
            }
        }
        x[i] = (y[i] - sigma) / diag; //updating x
    }
}

template <class T>
void CSRMatrix<T>::CholeskySolve(Matrix<T>* A, T* b, T* x)
{   
    // Check our dimensions match
    if (this->cols != this->rows)
    {
        std::cerr << "Ensure dimensions (rows and columns) of Matrix A match" << std::endl;
    }
    else {
        //Chelosky decomposing our Matrix A
        CSRMatrix<T>* Chol = this->CholeskyDecomp();

        //initialising y and filling it with 0s
        auto* y = new T[this->rows];

        for (int i = 0; i < this->rows; i++)
        {
            y[i] = 0;
        }

        //FORWARD SUBSTITUTION
        Chol->fsubstitution(b, y);


        //Transpose the Cholesky Matrix
        Chol->transposeiflower();

        //BACKWARD SUBSTITUTION
        Chol->bsubstitution(y, x);

        //No need to redeclare diag or sigma
    }
}