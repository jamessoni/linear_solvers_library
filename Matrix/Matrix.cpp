#include <iostream>
#include <math.h>
#include "Matrix.h"
#include <float.h> //For DBL_MAX
#include <numeric>
#include <vector>

using namespace std;

void daxpy(int n, double alpha, double *dx, int incx, double *dy, int incy)
{

   // Let's ignore the incx and incy for now
   /////#pragma loop(no_vector)
   for (int i = 0; i < n; i++)
   {
      dy[i * incy] += alpha * dx[i * incx];
   }
}

void daxpytx(int n, double alpha, double *dx, int incx, double *dy, int incy)
{
   //daxpy but the result is stored in x instead
   for (int i = 0; i < n; i++)
   {
      dx[i * incx] = dy[i*incy] + alpha * dx[i * incx];
   }
}

void dcopy(int n, double *dx, int incx, double *dy, int incy)
{
   for (int i = 0; i < n; i++)
   {
      dy[i * incy] = dx[i * incx];
   }

}

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
                this->values[i * cols + j] = (rand() % 2) + 0;
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
float Matrix<T>::RMS_norm_diff(T* vec_a, T* vec_b) 
{
    // RMS norm of the different of two vectors 
    // all input vectors/arrays need to be same size

    float sum_a = 0;

    // loop over all values in arraz
    for (int i = 0; i < this->rows; i++)
    {
        // add the squared difference to sum_a
        sum_a += (vec_a[i] - vec_b[i])*(vec_a[i] - vec_b[i]);

    }

    // return RMS norm of the squared difference
    return sqrt(sum_a / this->rows);
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
void Matrix<T>::gauss_seidel(Matrix<T>& a, Matrix<T>& b, Matrix<T>& x_init)
{
    //   Gauss-seidel implementation
    //   Method for solving a linear system, Ax = b, where A is a positive definite matrix
    //   Both convergence tolerance and fixed iteration methodologies presented 
    //   as convergence criteria
    //   Convergence tolerance and iteration number are predefined in the variables tol and iter_max below.

    double tol = 1e-5;
    int iter_max = 500;
    int iter = 0;
    double conve = 10;

    T* pout2 = new T[x_init.rows];

    for (int i = 0; i < x_init.rows * x_init.cols; i++)
    {
        x_init.values[i] = 0;
    }

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
            x_init.values[i] = b.values[i] / a.values[i * rows + i];
            for (int j = 0; j < this->cols; j++)
            {
                if (i == j)
                {
                    continue;
                }
                x_init.values[i] = x_init.values[i] - ((a.values[i * rows + j] / a.values[i * rows + i]) * x_init.values[j]);
            }
        }
        conve = RMS_norm_diff(pout2, x_init.values);
        std::cout << "\nconvergence values: " << conve;
    }
    std::cout << std::endl;

    delete[] pout2;

    //Same implementation as above using a predefined iteration number
    //while (iter_max > 0)
    //{
    //    
    //    for (int i = 0; i < this->rows; i++)
    //    {
    //        x_init.values[i] = b.values[i] / a.values[i * rows + i];
    //        for (int j = 0; j < this->cols; j++)
    //        {
    //            if (i == j)
    //            {
    //                continue;
    //            }
    //           x_init.values[i] = x_init.values[i] - ((a.values[i * rows + j] / a.values[i * rows + i]) * x_init.values[j]);
    //        }
    //    }   
    //    iter -= 1;
    //}  
    //std::cout << std::endl;

     //delete[] pout2;
     
}


template <class T>
void Matrix<T>::jacobi_solver_element(T* b, T* output, int maxIter, bool initialised) {
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
    int n = 0;
    double sum_RMS = 0;
    //set solution tolerance to e-10
    double tol = 1.e-10;


    // if not initialised fist input than use random numbers to 
    // initialise x_{k}
    if (initialised == false) {
        // initialise starting condition x_{k} 
        for (int i = 0; i < this->rows; i++)
        {
            output[i] = rand() % 50 + 5;
        }
    }

    // initialise starting condition x_{k} 
    for (int i = 0; i < this->rows; i++)
    {
        output[i] = rand() % 50 + 5;
    }

    // start iteration, only do maxIter steps 
    for (int n = 0; n < maxIter; n++) {
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

        // set RMS norm summation varible to zero
        sum_RMS = 0;

        // In this loop the RMS residual is calculated just like in the RMS_norm_diff func
        // but this loop is also used to copy the array therefore only one loop is used below
        // rather than two seperate loops
        for (int i = 0; i < this->rows; i++) {


            // add the squared difference to sum_a for the RMS convergence between the
            // different itterations
            sum_RMS += (pout2[i] - output[i]) * (pout2[i] - output[i]);

            // copy values into new array for next itteration
            output[i] = pout2[i];

        }

        // RMS norm of the squared difference
        conve = sqrt(sum_RMS / this->rows);

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
void  Matrix<T>::jacobi_solver_matrix(double* b, double* xk1, int maxIter, bool initialised) {
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

    auto* N = new Matrix<T>(this->rows, this->cols, true);
    auto* D = new Matrix<T>(this->rows, this->cols, true);



    T* xk2 = new T[this->rows];   // store x_{k+1}
    T* pvecvecarray = new T[this->rows];   // store vecVecsubtract
    T* pmatvecarray = new T[this->rows];   // store matVecMult
    //set solution tolerance to e-10
    double tol = 1.e-10;

    // initialize conve varible for first itteration 
    double conve = 12;

    // decompose value into a inverse diagonal D, and lower/upper triangular N
    this->jacobi_decomposition(D, N);

    // if not initialised input, use random numbers to 
    // initialise x_{k}
    if (initialised == false) {
        // initialise starting condition x_{k} 
        for (int i = 0; i < this->rows; i++)
        {
            xk1[i] = rand() % 50 + 5;
        }
    }

    // x_{k+1} = D(b - N * x_k)
    for (int n = 0; n < maxIter; n++) {
        N->matVecMult(xk1, pmatvecarray); // K = N * x_n
        N->vecVecsubtract(b, pmatvecarray, pvecvecarray); // R = b - K

        // need to copy old solution before it is overwritten in next step
        // this is needed to find RMS norm of the previous solution to the
        // current one
        for (int i = 0; i < this->rows; i++)
        {
            xk2[i] = xk1[i];
        }

        D->matVecMult(pvecvecarray, xk1); // X_{n+1} = D * R

        // find RMS norm of output-pout2 array 
        conve = RMS_norm_diff(xk2, xk1);

        // if rms norm is smaller than tolerance -> break loop
        if (conve < tol) {
            break;
        }
    }


   //delete[] N->values;
   //delete[] D->values;
   delete N;
   delete D;
   delete[] pvecvecarray;
   delete[] pmatvecarray;
   delete[] xk2;


}


template <class T>
void Matrix<T>::LUDecomp(Matrix<T>& L, Matrix<T>& U)
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
      double factor = 0; //factor by which to multiply the pivot row before adding it
      for(int i = step+1; i < U.rows; i++)
      {
         factor = U.values[(i)*U.cols+step]/U.values[step*U.cols+step];
         daxpy(U.cols-step,-factor,&U.values[step*U.cols+step],1,&U.values[i*U.cols+step],1);
         L.values[i*cols+step] = factor;
      }
   }
}

template <class T>
void Matrix<T>::SLUDecomp(Matrix<T>* LU)
{
   //Decompose to a single matrix

   // Check our dimensions match
   if (this->cols != LU->cols || this->cols != LU->cols)
   {
      std::cerr << "L and U must be of the same size as the matrix you want to decompose" << std::endl;
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
   for (int step = 0; step<LU->rows; step++)
   {
      double factor = 0;
      for(int i = step+1; i < LU->rows; i++)
      {
         factor = LU->values[(i)*LU->cols+step]/LU->values[step*LU->cols+step];
         daxpy(LU->cols-step,-factor,&LU->values[step*LU->cols+step],1,&LU->values[i*LU->cols+step],1);
         LU->values[i*cols+step] = factor;
      }
   }
}

template <class T>
void Matrix<T>::IPLUDecomp()
{
   //In place decomposition

   // Check if our LU matrix has had space allocated to it
   if (this->values == nullptr) 
   {
      std::cout << "ERROR. CANNOT DECOMPOSE NON-EXISTING MATRIX";
   }

   for (int step = 0; step<this->rows; step++)
   {
      double factor = 0;
      for(int i = step+1; i < this->rows; i++)
      {
         factor = this->values[(i)*this->cols+step]/this->values[step*this->cols+step];
         daxpy(this->cols-step,-factor,&this->values[step*this->cols+step],1,&this->values[i*this->cols+step],1);
         this->values[i*cols+step] = factor;
      }
   }
}

template <class T>
void Matrix<T>::fsubstitution(Matrix<T>& L, T* y, T* b)
{  
    for (int i = 0; i<L.cols;i++)
   {
      double sum =0;
      for (int j = 0; j<L.cols; j++)
      {
         sum+= L.values[i*L.cols + j]*y[j];
      }
      y[i]=(b[i]-sum)/1;
   }
}

template <class T>
void Matrix<T>::bsubstitution(Matrix<T>& U, T* x, T* y)
{
   for (int i = U.cols-1; i>=0;i--)
   {
      double sum = 0;
      for (int j = 0; j<U.cols; j++)
      {
         sum+= U.values[i*U.cols + j]*x[j];
      }

      x[i]=(y[i]-sum)/U.values[i*(U.cols+1)];
   }
}

template <class T>
void Matrix<T>::LUSolve(double* b, double* output, bool inplace)
{
   //pivoting

   Matrix<T> *LU;

   if (inplace) //inplace decomp
   {
      this->IPLUDecomp();
      LU = this;
   }
   else //not inplace decomp
   {
      auto *NewLU = new Matrix<T>(this->rows, this->cols, true);
      this->SLUDecomp(NewLU);
      LU = NewLU;
   }

   for(int i = 0; i<LU->cols;i++)
   {
      output[i] = 0;
   }

   //creating y vector for our intermediate step.
   T y[LU->cols];

   //forward substitution
   fsubstitution(*LU, y, b);

   //backward substitution
   bsubstitution(*LU,output,y);
   
   
   if (inplace) //delete LU if new space was allocated
   {
      delete LU;
   }
}

template <class T>
void Matrix<T>::conjugate_gradient(T* b, T* x, int maxIter, float tol)
{
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

   
   bnum = inner_product(r,r+this->cols,r,0); //finding inner product of r
	error = sqrt(bnum / this->cols); //finding error, using RMS

	while (count < maxIter && tol < error)
	{
      //CALCULATE ALPHA = (r*r)/(p(Ap))
		this->matVecMult(p, Ap); // Finding Ap
      anum = inner_product(r,r+this->cols,r,0.0); //r*r inner product
		aden = inner_product(p,p+this->cols,Ap,0.0); //finding (p(Ap))
		alpha = anum / aden; // calculating alpha
      
		// UPDATE X. x = x + ap
      daxpy(this->cols,alpha,p,1,x,1);

		// UPDATE R. R = R - alpha*(Ap)
      dcopy(this->cols, r, 1, r_old, 1); // storing r for later use in bden. r_old = r
      daxpy(this->cols,-alpha,Ap,1,r,1);
 
		//Calculate Beta
      bnum = inner_product(r,r+this->cols,r,0.0); //r*r inner product
		bden = inner_product(r_old,r_old+this->cols,r_old,0.0);
		beta = bnum / bden; // calculating beta


		//Update p = r + beta *p(old)
      daxpytx(this->cols,beta,p,1,r,1);

	   //Update error
		error = sqrt(bnum/this->cols); // RMS error 

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