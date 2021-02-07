#include <iostream>
#include <math.h>
#include "Matrix.h"
#include <float.h> //For DBL_MAX
#include <numeric>
#include <vector>

using namespace std;

// Constructor - using an initialisation list here
template <class T>
Matrix<T>::Matrix(int rows, int cols, bool preallocate) : rows(rows), cols(cols), size_of_values(rows* cols), preallocated(preallocate)
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
Matrix<T>::Matrix(int rows, int cols, T* values_ptr) : rows(rows), cols(cols), size_of_values(rows* cols), values(values_ptr)
{}

//Constructor - creating the SPD matrix
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
                // updating Matrix values on the diagonal - 
                // enabling Diagonal Dominance
                this->values[i * cols + j] = (rand() % diag_max) + diag_min;
            }
            else
            {
                // updating Matrix off-diagonal elements -
                // ensuring symmetry
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
    if (this->preallocated) {
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
    /*  This method returns true if the matrix is Symmetric
        and Diagonal Dominance occurs (the sum of the non diagonal entries of the row
        is less than the smallest diagonal entry)
    */

 

    int sum_nondiag{};
    int diag_val{};
    bool SPD{ false };

 

    //criteria 1. testing for symmetry in matrix
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            if (i == j)
            {
                continue;
            }
            else {
                // updating SPD bool condition as we loop over symmetric
                // elements within the Matrix
                if (this->values[i * rows + j] == this->values[j * cols + i]) {
                    SPD = true;
                }
                else {
                    SPD = false;
                }
            }
        }
    }
    //criteria 2. SPD test (weak) that enables convergence to occur 
    //            with certain solvers
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            if (i == j) {
                //ensuring Diagonal 
                diag_val = this->values[i * rows + j];
                continue;
            }
            else {
                sum_nondiag += this->values[i * rows + j];
            }
            if (diag_val <= sum_nondiag) {
                SPD = false;
            }
        }
    }
    if (SPD == true) {
        std::cout << "\nPasses weak SPD check";
    }
    else {
        std::cout << "\nMay not converge with chosen solvers";
    }
    return SPD;
}

template <class T>
void Matrix<T>::matVecMult(T* vec, T* output)
{
    // Matrix vector prodcut M * b = c
    // input vec is an array which is vector b
    // output in a array of size matrix->rows to store values
     
    // Set values to zero
    for (int i = 0; i < this->rows; i++)
    {
        output[i] = 0;
    }

    // loop over rows and cols
    for (int i = 0; i < this->rows; i++) 
    {
        for (int j = 0; j < this->cols; j++) 
        {
            // add the right values to the output array
            output[i] += this->values[i * this->cols + j] * vec[j];
        }
    }
}


template <class T>
void Matrix<T>::vecVecsubtract(T* vec_a, T* vec_b, T* output) 
{
    // Vector vector subtraction Vec_A - Vec_b = output
    // all three vectors/arrays need to be same size
  
    for (int i = 0; i < this->rows; i++) 
    {

        output[i] = vec_a[i] - vec_b[i];
    }
}

template <class T>
float Matrix<T>::RMS_norm_diff(T* vec_a, T* vec_b)
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
void Matrix<T>::gauss_seidel(Matrix<T>* A, T* b, T* x, int maxIter, float tol)
{

    /*
    Gauss-seidel solver implementation
        Solves a linear system of equations A * x = b using an iterative apporach
        Input :
            <T>array[] * a : SPD Matrix input
            <T>array[] * b : RHS of the linear system
            <T>array[] * x_init : array in which the solution will be stored in
            double tol : tolerance of the solver

        Output :
            none

        A needs to be a Diagonal Dominant / SPD matrix with no zeros on main diagonal
        and the linear system needs to have a solution.
        For gauss-seidel to be ran on maximum iterations - uncomment section below.
        int iter_max : maximum number of iterations can be altered below

    */
    // Check our dimensions match
    if (this->cols != this->rows)
    {
        std::cerr << "Ensure dimensions (rows and columns) of Matrix A match" << std::endl;
    }
    else {
        double RMS = 10;
        int count = 0;
        double* answer_check = new double[this->rows];
        T* pout2 = new T[this->rows];

        for (int i = 0; i < this->rows; i++)
        {
            x[i] = 0;
        }
        
        while (count < maxIter)
        {
            count += 1;
            //std::cout << "\niteration: " << iter;
            //updating pout2 as previous iteration x.values
            //enables convergence paramter to be checked against predefined tolerance
            for (int i = 0; i < this->rows; i++)
            {
                pout2[i] = x[i];
            }
            //looping of rows of A
            for (int i = 0; i < this->rows; i++)
            {
                x[i] = b[i] / this->values[i * this->rows + i];
                //looping over cols of A
                for (int j = 0; j < this->cols; j++)
                {
                    //skips over summation where (i==j)
                    if (i == j)
                    {
                        continue;
                    }
                    //computing x[i] = x[i] - (A[i][j] / A[i][i])*x[j]
                    x[i] = x[i] - ((A->values[i * this->rows + j] / A->values[i * this->rows + i]) * x[j]);
                }
            }
            // find error after every itteration
            A->matVecMult(x, answer_check);
            RMS = A->RMS_norm_diff(b, answer_check);
            if (RMS < tol) 
            {
                break;
            }
        }

        delete[] answer_check;
        delete[] pout2;
        /*
        // Gauss-seidel implementation with convergence criteria
        // using a predefined iteration number iter_max

        int iter_max = 500;
        while (iter_max > 0)
        {

            for (int i = 0; i < this->rows; i++)
            {
                x[i] = b[i] / a[i * rows + i];
                for (int j = 0; j < this->cols; j++)
                {
                    if (i == j)
                    {
                        continue;
                    }
                   x[i] = x[i] - ((a[i * rows + j] / a[i * rows + i]) * x[j]);
                }
            }
            iter -= 1;
        }

        delete[] pout2;
        */
    }
}


template <class T>
void Matrix<T>::jacobi_solver_element(Matrix<T>* A,T* b, T* x, int maxIter, bool initialised, float tol) {
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
    // Check our dimensions match
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
        int n = 0;
        double RMS = 0;
        double* answer_check = new double[A->rows];


        // if not initialised fist input than use random numbers to 
        // initialise x_{k}
        if (initialised == false) {
            // initialise starting condition x_{k} 
            for (int i = 0; i < A->rows; i++)
            {
                x[i] = rand() % 50 + 5;
            }
        }

        // initialise starting condition x_{k} 
        for (int i = 0; i < A->rows; i++)
        {
            x[i] = rand() % 50 + 5;
        }

        // start iteration, only do maxIter steps 
        for (int n = 0; n < maxIter; n++) {
            // loop over rows
            for (int i = 0; i < A->rows; i++) {
                // set variable to zero such that it can be added to 
                sum = 0;
                // loop over cols of matrix
                for (int j = 0; j < A->cols; j++) {
                    // if i = j dont do anything because that is row value that is calcualted
                    // for i not equal to j, mutiply both values and add to sum
                    if (i != j) {
                        sum += A->values[i * A->cols + j] * x[j];
                    }

                }
                // put the sum value in the right spot in the array with some additonal calculations
                pout2[i] = (1 / (A->values[i * A->rows + i])) * (b[i] - sum);
            }

            // loop is used to copy the new solution into the other array
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
void  Matrix<T>::jacobi_solver_matrix(Matrix<T>* A, T* b, T* x, int maxIter, bool initialised, float tol) {
    /*
    Jacobi solver using matrix manipulation
    Solves a linear system of equations A*x=b using an ittertive apporach, where every itteration
    mutiple matrix operations are performed: xk2 = D*(b - N * xk1)
    Input:
        <T>array[]* b: RHS of the linear system
        <T>array[]* xk1: array in which the solution will be stored in
        int maxIter: maxiumum itterations
        bool initialised: If true, output array was initialised with fist guess
                          If false, output array will be filled with random values
    Output:
        none

    A needs to be a SPD matrix wiht no zeros on main diagonal
    and the linear system needs to have a solution.
    Tolerance of the solver can be changed in the tol variable below

    */
    // Check our dimensions match
    // Check our dimensions match
    if (this->cols != this->rows)
    {
        std::cerr << "Ensure dimensions (rows and columns) of Matrix A match" << std::endl;
    }
    else {
        auto* N = new Matrix<T>(A->rows, A->cols, true);
        auto* D = new Matrix<T>(A->rows, A->cols, true);



        T* xk2 = new T[A->rows];   // store x_{k+1}
        T* pvecvecarray = new T[A->rows];   // store vecVecsubtract
        T* pmatvecarray = new T[A->rows];   // store matVecMult

        // initialize conve varible for first itteration 
        double RMS = 12;
        double* answer_check = new double[A->rows];

        // decompose value into a inverse diagonal D, and lower/upper triangular N
        A->jacobi_decomposition(D, N);

        // if not initialised input, use random numbers to 
        // initialise x_{k}
        if (initialised == false) {
            // initialise starting condition x_{k} 
            for (int i = 0; i < A->rows; i++)
            {
                x[i] = rand() % 50 + 5;
            }
        }

        // x_{k+1} = D(b - N * x_k)
        for (int n = 0; n < maxIter; n++) {
            N->matVecMult(x, pmatvecarray); // K = N * x_n
            N->vecVecsubtract(b, pmatvecarray, pvecvecarray); // R = b - K

            // need to copy old solution before it is overwritten in next step
            // this is needed to find RMS norm of the previous solution to the
            // current one
            for (int i = 0; i < A->rows; i++)
            {
                xk2[i] = x[i];
            }

            D->matVecMult(pvecvecarray, x); // X_{n+1} = D * R

          // find error after every itteration
            A->matVecMult(x, answer_check);
            RMS = A->RMS_norm_diff(b, answer_check);
            // if rms norm is smaller than tolerance -> break lool
            if (RMS < tol) {
                break;
            }
        }


        //delete[] N->values;
        //delete[] D->values;
        delete N;
        delete D;
        delete[] answer_check;
        delete[] pvecvecarray;
        delete[] pmatvecarray;
        delete[] xk2;
    }
}


template <class T>
void Matrix<T>::LUDecomp(Matrix<T>& L, Matrix<T>& U)
{
    //NOTE: UNSTABLE BECAUSE THERE IS NO PARTIAL PIVOTING

    // Check our dimensions match
    if (this->cols != L.cols || this->cols != U.cols || this->rows != L.rows || this->rows != U.rows )
    {
        std::cerr << "L and U must be square and of the same size as the matrix you want to decompose" << std::endl;
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
        for (int j = 0; j < this->rows; j++)
        {
            if (i == j)
            {
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

    //looping through every row
    for (int step = 0; step<U.rows; step++)
    {
        //factor by which to multiply the pivot row before adding it
        double factor = 0; 

        //Gaussian elimination, row by row (rows below the current pivot).
        for(int i = step+1; i < U.rows; i++)
        {
            //factor = value in pivot column of row being considered/diagonal value
            factor = U.values[(i)*U.cols+step]/U.values[step*U.cols+step];

            //daxpying  -factor*(pivot row) + (current row)
            daxpy(U.cols-step,-factor,&U.values[step*U.cols+step],1,&U.values[i*U.cols+step],1);

            //adding our factor to the L matrix
            L.values[i*cols+step] = factor;
        }
    }
}

template <class T>
void Matrix<T>::SLUDecomp(Matrix<T>* LU)
{
    //NOTE: UNSTABLE BECAUSE THERE IS NO PARTIAL PIVOTING

    // Check our dimensions match
    if (this->cols != LU->cols || this->rows != LU->rows)
    {
        std::cerr << "LU must be square and of the same size as the matrix you want to decompose" << std::endl;
        return;
    }

    // Check if our LU matrix has had space allocated to it
    if (LU->values != nullptr) 
    {
        LU->values = new T[this->rows * this->cols];
        // Don't forget to set preallocate to true now it is protected
        LU->preallocated = true;
    }

    // LU = A as a starting point
    for (int i=0;i<this->size_of_values; i++)
    {
        LU->values[i] = this->values[i];
    }

    //looping through every row
    for (int step = 0; step<LU->rows; step++)
    {
        //factor by which to multiply the pivot row before adding it
        double factor = 0; 

        //Gaussian elimination, row by row (rows below the current pivot).
        for(int i = step+1; i < LU->rows; i++)
        {
            //factor = value in pivot column of row being considered/diagonal value
            factor = LU->values[(i)*LU->cols+step]/LU->values[step*LU->cols+step];

            //daxpying  -factor*(pivot row) + (current row)
            daxpy(LU->cols-step,-factor,&LU->values[step*LU->cols+step],1,&LU->values[i*LU->cols+step],1);

            //adding our factor to the LU matrix
            LU->values[i*cols+step] = factor;
        }
    }
}

template <class T>
void Matrix<T>::fsubstitutionLU(Matrix<T>& L, T* y, T* b)
{  
    //Looping through every row
    for (int i = 0; i<L.cols;i++)
    {
        //Keep track of the total product of the row with the construction-in-progress y
        double sum = 0;

        //Loop through elements in the row
        for (int j = 0; j<L.cols; j++)
        {
            sum+= L.values[i*L.cols + j]*y[j];
        }

        //Update y, divide by 1 instead of the diagonal element because it's (theoretically) LU even though it's stored in a single matrix
        y[i]=(b[i]-sum)/1;
    }
}

template <class T>
void Matrix<T>::fsubstitution(Matrix<T>& L, T* y, T* b)
{  
    //Looping through every row
    for (int i = 0; i<L.cols;i++)
    {
        //Keep track of the total product of the row with the construction-in-progress y
        double sum = 0;

        //Loop through elements in the row
        for (int j = 0; j<L.cols; j++)
        {
            sum+= L.values[i*L.cols + j]*y[j];
        }

        //Update y, divide by the diagonal element
        y[i]=(b[i]-sum)/L.values[i*(L.cols+1)];
    }
}

template <class T>
void Matrix<T>::bsubstitution(Matrix<T>& U, T* x, T* y)
{
    //Looping through every row, backwards
   for (int i = U.cols-1; i>=0;i--)
   {
       //Keep track of the total product of the row with the construction-in-progress y
      double sum = 0;

      //Loop through elements in the row
      for (int j = 0; j<U.cols; j++)
      {
         sum+= U.values[i*U.cols + j]*x[j];
      }
    
    //Update x, divide by the diagonal element
      x[i]=(y[i]-sum)/U.values[i*(U.cols+1)];
   }
}

template <class T>
void Matrix<T>::LUSolve(Matrix<T>* A, T* b, T* x, bool inplace)
{
    // Check our dimensions match
    if (this->cols != this->rows)
    {
        std::cerr << "Ensure dimensions (rows and columns) of Matrix A match" << std::endl;
    }
    else {
        Matrix<T>* LU;

        if (inplace) //inplace decomp
        {
            this->SLUDecomp(this);
            LU = this;
        }
        else //not inplace decomp
        {
            auto* NewLU = new Matrix<T>(this->rows, this->cols, true);
            this->SLUDecomp(NewLU);
            LU = NewLU;
        }

        for (int i = 0; i < LU->cols; i++)
        {
            x[i] = 0;
        }

        //creating y vector for our intermediate step.
        //T y[LU->cols];
        auto* y = new T[this->cols];

        for (int i = 0; i < LU->cols; i++)
        {
            y[i] = 0;
        }

        //forward substitution
        fsubstitutionLU(*LU, y, b);

        //backward substitution
        bsubstitution(*LU, x, y);

        if (inplace) //delete LU if new space was allocated
        {
            delete LU;
        }
    }
}

template <class T>
void Matrix<T>::conjugate_gradient(Matrix<T>* A, T* b, T* x, int maxIter, float tol)
{
    // Check our dimensions match
    if (this->cols != this->rows)
    {
        std::cerr << "Ensure dimensions (rows and columns) of Matrix A match" << std::endl;
    }
    else {

        int count = 0; //iteration counter
        // PREALLOCATING EVERYTHING. DOES NOT MAKE SENSE TO ALLOCATE EACH LOOP.
        auto* r = new T[this->cols]; // residual r
        auto* Ax = new T[this->cols]; // array for storing Ax when calculating residual
        auto* p = new T[this->cols]; // direction p 
        auto* Ap = new T[this->cols]; // array for storing Ap when calculating alpha
        auto* r_old = new T[this->cols]; // r sized array to store the old r when calculating beta

        double alpha; //alpha
        double anum = 0; // for  r*r dot product, numerator for alpha coefficient 
        double aden = 0; // denominator for alpha coefficient

        double beta; //beta
        double bnum = 0; // for r*r dot product, numerator for beta coefficient
        double bden = 0; // for p(Ap) , denominator for beta coefficient

        double error = 0;  // error


        //// CALCULATING INITIAL r, P AND ERROR

       //calculate residual r = b - Ax
        this->matVecMult(x, Ax); // Find Ax
        this->vecVecsubtract(b, Ax, r); //r = b - Ax
        dcopy(this->cols, r, 1, p, 1); // first p = r


        bnum = inner_product(r, r + this->cols, r, 0); //finding inner product of r
        error = sqrt(bnum / this->cols); //finding error, using RMS

        while (count < maxIter && tol < error)
        {
            //CALCULATE ALPHA = (r*r)/(p(Ap))
            this->matVecMult(p, Ap); // Finding Ap
            anum = inner_product(r, r + this->cols, r, 0.0); //r*r inner product
            aden = inner_product(p, p + this->cols, Ap, 0.0); //finding (p(Ap))
            alpha = anum / aden; // calculating alpha

            // UPDATE X. x = x + ap
            daxpy(this->cols, alpha, p, 1, x, 1);

            // UPDATE R. R = R - alpha*(Ap)
            dcopy(this->cols, r, 1, r_old, 1); // storing r for later use in bden. r_old = r
            daxpy(this->cols, -alpha, Ap, 1, r, 1);

            //Calculate Beta
            bnum = inner_product(r, r + this->cols, r, 0.0); //r*r inner product
            bden = inner_product(r_old, r_old + this->cols, r_old, 0.0);
            beta = bnum / bden; // calculating beta


            //Update p = r + beta *p(old)
            daxpytx(this->cols, beta, p, 1, r, 1);

            //Update error
            error = sqrt(bnum / this->cols); // RMS error 

            // Add 1 to the iteration count
            count++;

            if (count == maxIter)
            {
                break;
            }
        }

        delete[] r;
        delete[] r_old;
        delete[] Ax;
        delete[] p;
        delete[] Ap;
    }
}

template<class T>
void Matrix<T>::CholeskyDecomp(Matrix<T>* L)
{  
    //REMEMBER TO ADD A CHECK FOR SAME NUMBER OF COLS AND ROWS
    // Check if our L matrix has had space allocated to it
    if (L->values != nullptr) 
    {
        L->values = new T[this->rows * this->cols];
        // Don't forget to set preallocate to true now it is protected
        L->preallocated = true;
    }

    //Decomp
    for (int i = 0; i < this->rows; i++)
    {  
        //all the non diagonal elements first
        for (int j = 0; j < i; j++)
        {
            double sigma = 0;
            // non-diagonals L(i, j) 
            for (int p = 0; p < j; p++)
            {
                sigma += L->values[i*this->cols + p] * L->values[j*this->cols + p]; 
            }
            L->values[i*this->cols + j] = ((this->values[i*this->cols + j] - sigma) / L->values[j*this->cols + j]);
        }

        //Then the diagonal element
        double sigma = 0;
        for (int p = 0; p < i; p++) 
        {
            sigma += L->values[i*this->cols + p] * L->values[i*this->cols + p];
        }

        L->values[i*this->cols + i] = sqrt(this->values[i*this->cols + i] - sigma);

        for (int j = i+1; j < this->cols; j++)
        {
            L->values[i*this->cols + j] = 0;
        }
    }
}

template<class T>
void Matrix<T>::transpose()
{
   //slot we can store one of our values in while we switch them
   T keepval;

   //looping over ij, Aij=Aji
   for (int i = 0; i < this->rows; i++)
   {
      for (int j = 0; j < i; j++)
      {
         keepval = this->values[i*this->cols+j];
         this->values[i*this->cols+j] = this->values[j*this->cols+i];
         this->values[j*this->cols+i] = keepval;
      }
   }
}

template<class T>
void Matrix<T>::CholeskySolve(Matrix<T>* A, T* b, T* x)
{  
    // Check our dimensions match
    if (this->cols != this->rows)
    {
        std::cerr << "Ensure dimensions (rows and columns) of Matrix A match" << std::endl;
    }
    else {
        //Decomp
        auto* L = new Matrix<double>(this->rows, this->cols, true);
        this->CholeskyDecomp(L);

        //creating an intermediate y;
        //T y[this->cols];
        T* y = new double[this->cols];

        //filling y with 0s
        for (int i = 0; i < this->rows; i++)
        {
            y[i] = 0;
        }


        //forward substitution
        fsubstitution(*L, y, b);

        //transposing our L
        L->transpose();

        //backward substitution
        bsubstitution(*L, x, y);

        delete L;
    }
}

template <class T>
void Matrix<T>::daxpy(int n, double alpha, double* dx, int incx, double* dy, int incy)
{

    // Let's ignore the incx and incy for now
    /////#pragma loop(no_vector)
    for (int i = 0; i < n; i++)
    {
        dy[i * incy] += alpha * dx[i * incx];
    }
}

template <class T>
void Matrix<T>::daxpytx(int n, double alpha, double* dx, int incx, double* dy, int incy)
{
    //daxpy but the result is stored in x instead
    for (int i = 0; i < n; i++)
    {
        dx[i * incx] = dy[i * incy] + alpha * dx[i * incx];
    }
}

template <class T>
void Matrix<T>::dcopy(int n, double* dx, int incx, double* dy, int incy)
{
    for (int i = 0; i < n; i++)
    {
        dy[i * incy] = dx[i * incx];
    }

}

template <class T>
void Matrix<T>::restrictions(Matrix<T>* I) {
   // create restriction matrix for the multigrid method

   //loop over all rows and cols and fill with zero
    for (int i = 0; i < I->rows; i++) {
        for (int j = 0; j < I->cols; j++) {
            
                I->values[i * I->cols + j] = 0;
            
        }
    }
    //loop over all rows and cols and fill diagonals with 0.5 and the two
    //off diagonals with 0.25
    for (int i = 0; i < I->rows; i++) {
        for (int j = 0; j < I->cols; j++) {
            if (i == j) {
                I->values[i * I->cols + j] = 0.25;
                I->values[i * I->cols + j + 1] = 0.5;
                I->values[i * I->cols + j + 2] = 0.25;

            }

        }
    }

    I->values[I->rows - 1] = 0;
    I->values[I->cols*I->rows - I->cols] = 0;
}

template <class T>
void Matrix<T>::prolongation(Matrix<T>* a, Matrix<T>* out) {
   //create prolongation matirx by transposing it

    for (int i = 0; i < a->rows; i++) {
        for (int j = 0; j < a->cols; j++) {
            (out->values[j * out->cols + i]) = 2*(a->values[i * a->cols + j]);
        }
    }

}



template <class T>
void Matrix<T>::vecVecadd(T* vec_a, T* vec_b) {
    // Vector vector addition Vec_A + Vec_b = output
    // all three vectors/arrays need to be same size

    for (int i = 0; i < this->rows; i++) {

        vec_a[i] +=  vec_b[i];

    }

}



template <class T>
void  Matrix<T>::onestep_multigrid(Matrix<T>* a, double* b, double* x, int maxIter, bool initialised)
{    /*
     Multigrid method to solve the 1-D Poisson equation: -\nabla^2 u = g using a central difference scheme in space.
     This will lead to the following stencil: (-u_{i+1} + 2u_i - u_{i-1})/(\deta x^2) which will be expressed in a matrix form:
            [2  -1  0  0  0 ] [u_{0}]   [g_{0}]
            [-1  2  -1  0  0] [u_{1}]   [g_{1}]
            [0  -1  2  -1  0] [u_{2}] = [g_{2}]
            [0  0  -1  2  -1] [u_{3}]   [g_{3}]
            [0  0  0  -1   2] [u_{4}]   [g_{4}]
    This linear system can be sovled using the general appraches i.e. jacobi or LU solvers but a multigrid method will be more effective. 
    
    Source: Murray P, Solving Linear ODEs With Multigrid [Internet]. 
    Available from: https://peytondmurray.github.io/coding/multigrid-solve-odes/#moving-between-grids
    
    */

    //declare varaibles
    auto* restrict_matrix = new Matrix<T>((a->rows - 1)/2, a->cols, true); // restrictor which is  (n-1)/2 x n
    auto* prolongation_matrix = new Matrix<T>(a->rows, (a->cols - 1)/2, true); // prolonation matrix which is n x (n-1)/2
    // matrices below are just for storing different vectors and matrices operations
    auto* A_P = new Matrix<T>(a->rows, (a->cols - 1) /2, true);
    auto* Ac = new Matrix<T>((a->cols - 1)/ 2, (a->cols - 1) / 2, true);
    double* Ax = new double[a->rows];
    double* r = new double[a->rows];
    double* b_c = new double[a->rows];
    

    // First, do smoothing steps on fine grid to reduce high frequency errors
    a->jacobi_solver_element(a, b, x, 8, false);

    // compute residual of fine grid r (r = b - A * x)
    a->matVecMult(x, Ax);
    a->vecVecsubtract(b, Ax, r);

    // create restrictor and prolongation matrix
    a->restrictions(restrict_matrix);
    restrict_matrix->values[1] = 0.5;
    prolongation_matrix->prolongation(restrict_matrix, prolongation_matrix);


    //do restriction Ac = restriction * A * prologation in two steps:
    // 1. A_P = A * prologation
    a->matMatMult(*prolongation_matrix, *A_P);

    // 2. Ac = restriction * A_P
    restrict_matrix->matMatMult(*A_P, *Ac);

   
    //restrict error to coarse mesh 
    restrict_matrix->matVecMult(r, b_c);

    // array for storing the solution of coarse smoothing 
    double* x_c = new double[Ac->rows];


     //smoothing on coarse gird until tolerance is reached
    a->jacobi_solver_matrix(Ac, b_c, x_c, 5000, false);


    // compute residual of coarse grid r (r = b_c - Ac * x_c)
    prolongation_matrix->matVecMult(x_c, r);

   
    // interpolating coarser grid error into fine grid.
    a->vecVecadd(x, r);


    // smoothing on fine grid 
    a->jacobi_solver_element(a, b, x, 8, true);


    // I know this code is leaking memory but it crashes when I use delete
    // have to use smart pointers next time 
    //delete smoothing_P;
    //delete smoothing_P_T;
    //delete A_P;
    //delete Ac;
    //delete[] Ax;
    //delete[] r;
    //delete[] r_c;

}