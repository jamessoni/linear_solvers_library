#include <iostream>
#include "Matrix.h"
#include <float.h> //For DBL_MAX
#include <numeric>

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
         //std::cout << "x: " << U.values[step*U.cols] << " y: "<< U.values[i*U.cols] << " n: " << U.cols-step << " factor: " << factor << "\n";
         
         daxpy(U.cols-step,-factor,&U.values[step*U.cols+step],1,&U.values[i*U.cols+step],1);
         
         //U.printMatrix();
         // for(int j = step; j < U.cols; j++)
         // {
         //    //U.values[i*U.cols+j] -= U.values[step*U.cols+j]*factor;
         // }         
         L.values[i*cols+step] = factor;
      }
   }
}

template <class T>
void Matrix<T>::SLUDecomp(Matrix& LU)
{
   //Decompose to a single matrix

   // Check our dimensions match
   if (this->cols != LU.cols || this->cols != LU.cols)
   {
      std::cerr << "L and U must be of the same size as the matrix you want to decompose" << std::endl;
      return;
   }

   // Check if our LU matrix has had space allocated to it
   if (LU.values != nullptr) 
   {
      LU.values = new T[this->rows * this->cols];
      // Don't forget to set preallocate to true now it is protected
      LU.preallocated = true;
   }

   // LU = A as a starting point
   for (int i=0;i<this->size_of_values; i++)
   {
      LU.values[i] = this->values[i];
   }

   for (int step = 0; step<LU.rows; step++)
   {
      double factor = 0;
      for(int i = step+1; i < LU.rows; i++)
      {
         factor = LU.values[(i)*LU.cols+step]/LU.values[step*LU.cols+step];
         daxpy(LU.cols-step,-factor,&LU.values[step*LU.cols+step],1,&LU.values[i*LU.cols+step],1);
         LU.values[i*cols+step] = factor;
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
      this->SLUDecomp(*NewLU);
      LU = NewLU;
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
		alpha = bnum / aden; // calculating alpha

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