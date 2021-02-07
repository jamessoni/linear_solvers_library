# Team Void ACSE-5 Group Assignment

### Implementing Linear Solvers

Within this Advanced Programming Assignment the project aimed to implement an algorithm to solve the linear system *A**x* = *b*, where A is a positive definite matrix, with x and b being vectors.
This linear solver library is written in the programming language, C++, maximising its advantages of Object Oriented Design, memory allocation and optimisation.
The assignment further evaluates, tests and aims to optimise the linear solvers created fulfilling convergence to solutions for both dense and sparse input matrices. Present within this linear solvers library, exists dense and sparse libraries - each with their respective .cpp and .h files. 

#### Getting Started / Pre-requisites
##### Mac/Linux:
IDE or C++ compiler of choice. Recommendations listed below:
* default g++ compiler
* Xcode
##### Windows users:
IDE or C++ compiler of choice. Recommendations listed below:
* Microsoft Visual Studio Community
* Visual Studio Code

##### Installation
Clone or download the [Team Void github repository](https://github.com/acse-2020/group-project-team-void) to your local machine

##### Using
For use of the program will depend on your OS and IDE/compiler of choice:
Visual Studio Community IDE users:
* Open a Visual Studio Project, add all the cloned files to the source files list
* Build project and run.

Terminal:
* Enter directory with cloned repository files within the Matrix/ folder

        $> g++ main_matrix.cpp -fpermissive
        $> .\a


#### Linear Solvers:
Within this library there exists various linear solvers; each with their enabled callable method within classes Matrix or CSRMatrix. They can all be accessed and used with the above main_matrix.cpp providing a UI platform to see examples as well as called specifically with the required inputs.

##### Implemented Solvers:
* Jacobi element solver
* Jacobi matrix solver (dense & sparse)
* Gauss-seidel solver (dense & sparse)
* LU solver
* Conjugate gradient solver
* Cholesky solver

#### Declaring input Matrices:
##### Dense Matrices:
The values of the dense matrices will depend on whether you as the user wish to input a specified matrix.
One example shown below:
'''
    int rows = 4;
    int cols = 4;

    //creating matrix
    auto *dense_mat = new Matrix<double>(rows, cols, true);

    //allocating specific matrix entries
    vector<double> vs = {5, 6, 2, 1, 9, 9, 7, 2, 4, 3, 8, 1, 2, 0, 9, 1};

    //filling Matrix dense_mat with values specified above
    for (int i = 0; i < rows * cols; i++)
    {
      dense_mat->values[i] = vs[i];
    }

    [...]

    delete dense_mat; 
'''
The example above declares a dense matrix which has further been filled with the values present in the vector (vs). 
However, in providing ease of use of the system - there is functuality to create a random SPD matrix if necessary.
Shown below:
'''
    int rows = 4;
    int cols = 4;

    //creating SPD matrix
    auto* A = new Matrix<double>(rows, cols, 3 * pow(rows, 2), 2 * pow(rows, 2));

    [...]

    delete A;
'''
This above with specific argument inputs calls class constructor; in turn producing a randomised Diagonal Dominant Positive Definite Matrix. This constructor call further ensures the matrix produced will converge with specified linear solvers adhering to both diagonal dominance and/or positive definiteness.

##### Sparse Matrices:
Further, sparse matrices like dense share familiar attributes and therefore call from both Matrix and sub-class CSRMatrix.
The sparse matrices are defined through compressed sparse row (CSR) format - pertaining to three One Dimensional arrays that respectively enables fast access of non zero values for methods and calculations.
Example of CSR below:
'''

    int rows = 4;
    int cols = 4;
    int const nnzs = 8;
    
    auto* sparse_mat = new CSRMatrix<double>(rows, cols, nnzs, true);

    int vals[nnzs] = { 10,2,12,2,6,1,1,9 };

    //filling sparse_mat values array
    for (int i = 0; i < nnzs; i++)
    {
        sparse_mat->values[i] = vals[i];
    }

    //filling sparse_mat col_index array
    int col_ind[nnzs] = { 0,2,1,0,2,3,2,3 };
    for (int i = 0; i < nnzs; i++)
    {
        sparse_mat->col_index[i] = col_ind[i];
    }

    //filling sparse_mat row_position array
    int row_ind[5] = { 0,2,3,6,8 };
    for (int i = 0; i < 5; i++)
    {
        sparse_mat->row_position[i] = row_ind[i];
    }

    [...]
    delete sparse_mat;
'''
##### Solvers
To call the linear solvers (for both dense and sparse) it will require using an object matrix. This will include calling a constructor which randomises its Matrix element values (e.g. SPD constructor) or a constructor where you wish to fill it with values.
Fulfilling the inputs of the solver (refer to the declared solver in header file), you can call a solver as follows:
###### creating inputs SPD Matrix, 1D arrays x and b:
'''

        float tol = 1e-6;
        int rows = 4;
        int cols = 4;
        
        //creating the SPD matrix
        auto* A = new Matrix<double>(rows, cols, 3 * pow(rows, 1.5), 2 * pow(rows, 1.5));

        //creating and filling array b with random values
        double* b = new double[rows];
        for (int i = 0; i < rows; i++) {
            b[i] = rand() % 300 + 5;
        }

        // set up and initialising x (solution) array
        double* x = new double[rows];
        for (int i = 0; i < rows; i++)
        {
            x[i] = 0.0;
        }
'''
###### Example of solver call:
'''
            
        //using gauss-seidel method
        A->gauss_seidel(*A, b, x, tol);
'''

###### Documentation:
The class Matrix has three constructors and a destructor:

'''
Matrix(int rows, int cols, bool preallocate);
Matrix(int rows, int cols, T *values_ptr);
Matrix(int rows, int cols, int diag_max, int diag_min);

virtual ~Matrix();
'''

The first constructor takes two ints for the number of rows and the number of columns and a further bool for the option to preallocate memory.
The second constructor takes two ints for the number of rows and the number of columns and a pointer to an array of values.

The third constructor is useful for testing: it takes two ints for the number of rows and the number of columns and two further integers diag_max and diag_min.
It creates a matrix with diagonal values randomised between diag_max and diag_min and all other entries = 1.

Example call:

'''
auto *dense_mat = new Matrix<double>(5, 5, 10,-10);
dense_mat->printMatrix();
-------------
Printing matrix

-9 1 1 1 1
1 -3 1 1 1 
1 1 -6 1 1
1 1 1 -10 1
1 1 1 1 -1
'''

Two methods to display the Matrix's values:

void printValues();
virtual void printMatrix();

A method for matrix-matrix multiplication. A.B=AB
It takes two matrices as input: one is the matrix to multiply (B) by and the other is the output matrix AB.
void matMatMult(Matrix<T>& mat_right, Matrix<T>& output);

Example:
'''
auto *dense_mat_left = new Matrix<double>(3, 3, 5,-5);
auto *dense_mat_right = new Matrix<double>(3, 3, 5,-5);
auto *dense_mat_output = new Matrix<double>(3, 3, true);
dense_mat_left->matMatMult(dense_mat_right,dense_mat_output);
-------------
Printing matrix

-4 1 1
1 -3 1
1 1 -1 
Printing matrix

-5 1 1
1 -1 1
1 1 -1
Printing matrix

22 -4 -4
-7 5 -3
-5 -1 3
'''

A method to check whether a Matrix is Symmetric Positive Definite.
The 
It returns a bool as well as printing a statement.
bool SPDMatrixcheck();

Example:
'''
auto *SPDM = new Matrix<double>(3, 3, 10,5);
SPDM->printMatrix();
bool spd = SPDM->SPDMatrixcheck();
std::cout << "\n"<< spd<< "\n";
SPDM->values[4] = 0;
SPDM->printMatrix();
spd = SPDM->SPDMatrixcheck();
std::cout << "\n" << spd<< "\n";
-------------
Printing matrix

6 1 1
1 12 1
1 1 9

Passes weak SPD check
1
Printing matrix

6 1 1
1 0 1
1 1 9

May not converge with chosen solvers
0
'''

A method for matrix vector multiplication Ax = b (in this case x is vec, b is output)
void matVecMult(T* vec, T* output);

Example:
'''
auto *Mat = new Matrix<double>(3, 3, true);

//Filling our matrix with values
double values[9] = {3,1,1,1,3,1,1,1,2};
for(int i = 0;i<9; i++)
{
    Mat->values[i]=values[i];
}

Mat->printMatrix();
double b[3] = {2,3,4};
double x[3];
Mat->matVecMult(b,x);

//printing x
for (int i = 0; i<3; i++)
{
    cout << x[i] << endl;
}
-------------
Printing matrix

3 1 1
1 3 1
1 1 2
13
15
13
'''

Two methods to perform vector vector subtraction and to calculate the RMS difference between two vectors.
void vecVecsubtract(T* vec_a, T* vec_b, T* output);
float RMS_norm_diff(T* vec_a, T* vec_b);

A few iterative solver methods for the system Ax = b. They all take A, b and x as inputs and also take an int maxIter and a float tol parameters.
The solver will continue iterating until it reaches maxIter iterations or until when the RMS difference between Ax and b is lower than tol.

The jacobi solver has a further bool: if it is set true, then it takes the x vector provided as a first guess. Otherwise it initialises it again.

void gauss_seidel(Matrix<T>* A, T* b, T* x, float tol);
void jacobi_solver_element(Matrix<T>* A,T* b, T* output, int maxIter, bool initialised, float tol);
void jacobi_solver_matrix(Matrix<T>* A, double* b, double* output, int maxIter, bool initialised, float tol);
void conjugate_gradient(Matrix<T>* A, T* b, T* x, int maxIter, float tol);

A few non iterative solver methods for the system Ax = b are also provided.
These generally achieve higher precisions, but take longer to execute.

They take A, b and x as inputs. The LUSolve method also 
void LUSolve(Matrix<T>* A, double* b, double* output, bool inplace);
void CholeskySolve(Matrix<T>* A, T* b, T* x);

void jacobi_decomposition(Matrix<T>* D, Matrix<T>* N);
void LUDecomp(Matrix<T>& L, Matrix<T>& U);
void SLUDecomp(Matrix<T>* LU);
void IPLUDecomp();
void fsubstitution(Matrix<T>& L, T* y,T* b);
void fsubstitutionLU(Matrix<T>& L, T* y,T* b);
void bsubstitution(Matrix<T>& U, T* x, T* y);
void transpose();



void CholeskyDecomp(Matrix<T>* L);

A few methods to "vectorise" operations in our other methods. They are included as Matrix methods and not as BLAS calls because it was giving us compilation errors
void daxpy(int n, double alpha, double* dx, int incx, double* dy, int incy);
void daxpytx(int n, double alpha, double* dx, int incx, double* dy, int incy);
void dcopy(int n, double* dx, int incx, double* dy, int incy);

The class CSRMatrix has the following methods:

CSRMatrix(int rows, int cols, int nnzs, bool preallocate);
CSRMatrix(int rows, int cols, int nnzs, T* values_ptr, int* row_position, int* col_index);
~CSRMatrix();

virtual void printMatrix();

void matVecMult(T* input, T* output);
void dense2sparse(Matrix<T>& tosparsify, CSRMatrix<T>* output);

void jacobi_solver_sparse(CSRMatrix<T>* A, T* b, T* output, int maxIter, bool initialised, float tol);

float RMS_norm_diff(T* vec_a, T* vec_b);

void gauss_seidel_sparse(CSRMatrix<T>* A, T* b, T* x_init, float tol);

int getv(int row,int col);
CSRMatrix<T>* matMatMult(CSRMatrix<T>& mat_right);

void transposeiflower();
CSRMatrix<T>* CholeskyDecomp();
void fsubstitution(T* b, T* y);
void bsubstitution(T* y, T* x);
void CholeskySolve(Matrix<T>* A, T* b, T* x);
        
#### License
See LICENSE.txt for details.

###### Team-Void ACSE-5
