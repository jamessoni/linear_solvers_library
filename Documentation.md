###### Documentation:
The class Matrix has three constructors and a destructor:

        Matrix(int rows, int cols, bool preallocate);
        Matrix(int rows, int cols, T *values_ptr);
        Matrix(int rows, int cols, int diag_max, int diag_min);

        virtual ~Matrix();

The first constructor takes two ints for the number of rows and the number of columns and a further bool for the option to preallocate memory.
The second constructor takes two ints for the number of rows and the number of columns and a pointer to an array of values.

The third constructor is useful for testing: it takes two ints for the number of rows and the number of columns and two further integers diag_max and diag_min.
It creates a matrix with diagonal values randomised between diag_max and diag_min and all other entries = 1.

Example call:

        auto *dense_mat = new Matrix<double>(5, 5, 10,-10);
        dense_mat->printMatrix();
        -------------
        Printing matrix

        -9 1 1 1 1
        1 -3 1 1 1 
        1 1 -6 1 1
        1 1 1 -10 1
        1 1 1 1 -1

Two methods to display the Matrix's values:

void printValues();
virtual void printMatrix();

A method for matrix-matrix multiplication. A.B=AB
It takes two matrices as input: one is the matrix to multiply (B) by and the other is the output matrix AB.
void matMatMult(Matrix<T>& mat_right, Matrix<T>& output);

Example:

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

A method to check whether a Matrix is Symmetric Positive Definite.
It returns a bool as well as printing a statement.
bool SPDMatrixcheck();

Example:

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

A method for matrix vector multiplication Ax = b (in this case x is vec, b is output)
void matVecMult(T* vec, T* output);

Example:

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

Two methods to perform vector vector subtraction and to calculate the RMS difference between two vectors.
void vecVecsubtract(T* vec_a, T* vec_b, T* output);

Example:

        auto* Mat = new Matrix<double>(4, 4, true);
        double vec_a[4]{ 10,2,4,9 };
        double vec_b[4]{ 2, 7,2,1 };
        double output[4]{};
        Mat->vecVecsubtract(vec_a, vec_b, output);

        //printing output
        for (int i = 0; i < 4; i++)
        {
            cout << output[i] << " ";
        }
        -------------
        Printing vector

        8 -5 2 8

float RMS_norm_diff(T* vec_a, T* vec_b);

Example:

        auto* A = new Matrix<double>(4, 4, 3 * pow(4, 2), 2 * pow(4, 2));

        // set up arrays to store solutions and answer check
        double* x = new double[4];
        double* b = new double[4];
        double* answer_check = new double[4];
        for (int i = 0; i < 4; i++) {
            b[i] = rand() % 300 + 5;
        }
        for (int i = 0; i < 4; i++)
        {
            x[i] = 0.0;
        }
        A->jacobi_solver_matrix(A, b, x, 1000, true, 1e-6);

        //if this is close to zero, the function works 
        double RMS = A->RMS_norm_diff(b, x);
        cout << "RMS: " << RMS;
        -------------
        Printing

        200.939

The transpose method takes in an object of the class matrix and performs Matrix transpose. A -> A<sup>T</sup>
void transpose();

Example:

        //creating Matrix
        auto* A = new Matrix<double>(4, 4, true);
        //filling the values of matrix A
        for (int i = 0; i < 4 * 4; i++)
        {
            A->values[i] = i;
        }
        A->printMatrix();
        A->transpose();
        A->printMatrix();
        -------------
        Printing matrix

        0 1 2 3
        4 5 6 7
        8 9 10 11
        12 13 14 15
        Printing matrix

        0 4 8 12
        1 5 9 13
        2 6 10 14
        3 7 11 15

The jacobi decomposition method will decompose A into two matrix D and N.
void jacobi_decomposition(Matrix<T>* D, Matrix<T>* N);

Example:

        //calling SPD function for random Matrix
        auto* A = new Matrix<double>(4, 4, 3 * pow(4, 2), 2 * pow(4, 2));
        //initialising matrices N and D
        auto* N = new Matrix<double>(A->rows, A->cols, true);
        auto* D = new Matrix<double>(A->rows, A->cols, true);
        A->jacobi_decomposition(D, N);
        A->printMatrix();
        D->printMatrix();
        N->printMatrix();
        -------------
        Printing matrix

        73 1 1 1
        1 67 1 1
        1 1 78 1
        1 1 1 36
        Printing matrix

        0.0136986 0 0 0
        0 0.0149254 0 0
        0 0 0.0128205 0
        0 0 0 0.0277778
        Printing matrix

        0 1 1 1
        1 0 1 1
        1 1 0 1
        1 1 1 0

These two methods decompose the input matrix into L and U, where SLUDecomp() stores L and U in a single matrix LU.
Unfortunately row pivoting has not been implemented and so the user needs to be careful.
void LUDecomp(Matrix<T>& L, Matrix<T>& U);
void SLUDecomp(Matrix<T>* LU);

Example:

        auto* A = new Matrix<double>(4, 4, 3 * pow(4, 2), 2 * pow(4, 2));
        auto* LU = new Matrix<double>(A->rows, A->cols, true);
        A->SLUDecomp(LU);
        A->printMatrix();
        LU->printMatrix();
        -------------
        Printing matrix

        73 1 1 1
        1 67 1 1
        1 1 78 1
        1 1 1 36
        Printing matrix

        73 1 1 1
        0.0136986 66.9863 0.986301 0.986301
        0.0136986 0.0147239 77.9718 0.971779
        0.0136986 0.0147239 0.0124632 35.9597

Cholesky decomposition decomposes the positive definite matrix into a lower triangular and its conjugate component.
void CholeskyDecomp(Matrix<T>* L);

Example:

        auto* A = new Matrix<double>(4, 4, 3 * pow(4, 2), 2 * pow(4, 2));
        auto* L = new Matrix<double>(A->rows, A->cols, true);
        A->CholeskyDecomp(L);
        A->printMatrix();
        L->printMatrix();
        -------------
        Printing matrix

        73 1 1 1
        1 67 1 1
        1 1 78 1
        1 1 1 36
        Printing matrix

        8.544 0 0 0
        0.117041 8.18452 0 0
        0.117041 0.120508 8.83016 0
        0.117041 0.120508 0.110052 5.99664

The three substitution method are used to solve the system LUx = b.
They assume an input of a lower and upper triangular matrix, respectively. Values above or below the diagonal will be ignored.
The fsubstitutionLU method assumes that diagonal elements are equal to 1 (as we are storing L and U on the same matrix).

void fsubstitution(Matrix<T>& L, T* y,T* b);
void fsubstitutionLU(Matrix<T>& L, T* y,T* b);
void bsubstitution(Matrix<T>& U, T* x, T* y);

A few methods to "vectorise" operations in our other methods. They are included as Matrix methods and not as BLAS calls because it was giving us compilation errors
void daxpy(int n, double alpha, double* dx, int incx, double* dy, int incy);
void daxpytx(int n, double alpha, double* dx, int incx, double* dy, int incy);
void dcopy(int n, double* dx, int incx, double* dy, int incy);

A few iterative solver methods for the system Ax = b. They all take A, b and x as inputs and also take an int maxIter and a float tol parameters.
The solver will continue iterating until it reaches maxIter iterations or until when the RMS difference between Ax and b is lower than tol.

The jacobi solver has a further bool: if it is set true, then it takes the x vector provided as a first guess. Otherwise it initialises it again.

void gauss_seidel(Matrix<T>* A, T* b, T* x, float tol);
void jacobi_solver_element(Matrix<T>* A,T* b, T* output, int maxIter, bool initialised, float tol);
void jacobi_solver_matrix(Matrix<T>* A, double* b, double* output, int maxIter, bool initialised, float tol);
void conjugate_gradient(Matrix<T>* A, T* b, T* x, int maxIter, float tol);

A few non iterative solver methods for the system Ax = b are also provided.
These generally achieve higher precisions, but take longer to execute.

LUsolve uses LU decomposition, Cholesky uses LL<sup>T</sup> decomposition, both using back/forward substitution methods.
The bool parameter inplace in LUSolve is to decide whether we want to allocate new memory to store the LU decomposed matrix or are happy to "destroy" our original matrix

void LUSolve(Matrix<T>* A, double* b, double* output, bool inplace);
void CholeskySolve(Matrix<T>* A, T* b, T* x);

The class CSRMatrix has the following methods:

These include two constructors which create an object of the class CSRMatrix.
The destructor will be called and automatically de-allocate memory itself if preallocated == true.
CSRMatrix(int rows, int cols, int nnzs, bool preallocate);
CSRMatrix(int rows, int cols, int nnzs, T* values_ptr, int* row_position, int* col_index);
~CSRMatrix();

virtual void printMatrix();

Method computing Sparse matrices with vectors product
void matVecMult(T* input, T* output);

Example:
    
    auto *sparse_mat = new CSRMatrix<double>(rows, cols, nnzs, true);
    // filling with i on diagonals
    for (int i = 0; i < nnzs; i++)
    {
    sparse_mat->values[i] = i;
    sparse_mat->col_index[i] = i;
    sparse_mat->row_position[i] = i;
    }
    // Now let's print it
    sparse_mat->printMatrix();
    // Now let's test our matvec
    double *input = new double[cols];
    double *output = new double[rows];
    for (int i = 0; i < cols; i++)
    {
        input[i] = i;
    }
    // computing matrix-vector product
    sparse_mat->matVecMult(input, output);
    // printing
    for (int i = 0; i < rows; i++)
    {
        cout << " " << output[i];
    }
    -------------
    Printing matrix
    Values: 0 1 2 3
    row_position: 0 1 2 3 -842150451
    col_index: 0 1 2 3
    0 1 4 0

Dense2sparse method takes in a dense matrix and "sparsifies" it.
void dense2sparse(Matrix<T>& tosparsify, CSRMatrix<T>* output);

Example:    

    int nnzs = 0;
    int rows = 3;
    int cols = 3;
    auto* A = new Matrix<double>(rows, cols, 3 * pow(rows, 2), 2 * pow(rows, 2));
    for (int j = 0; j < A->rows * A->cols; j++)
    {
        int val = A->values[j];
        if (val != 0)
        {
            nnzs += 1;
        }
    }
    //auto* sparse_mat = new CSRMatrix<double>(rows, cols, nnzs, true);
    A->printMatrix();
    auto* sparse_mat_2 = new CSRMatrix<double>(rows, cols, nnzs, true);
    sparse_mat_2->dense2sparse(*A, sparse_mat_2);
    sparse_mat_2->printMatrix();
    -------------
    Printing matrix

    32 1 1
    1 44 1
    1 1 34
    Printing matrix
    Values: 32 1 1 1 44 1 1 1 34
    row_position: 0 3 6 9
    col_index: 0 1 2 0 1 2 0 1 2

Two iterative solver methods for the system Ax = b, where A is a sparse matrix. They all take A, b and x as inputs and also take an int maxIter and a float tol parameters.
The solver will continue iterating until it reaches maxIter iterations or until when the RMS difference between Ax and b is lower than tol.
void jacobi_solver_sparse(CSRMatrix<T>* A, T* b, T* output, int maxIter, bool initialised, float tol);
void gauss_seidel_sparse(CSRMatrix<T>* A, T* b, T* x_init, float tol);

A method for matrix-matrix multiplication. A.B=AB
It takes two matrices as input: one is the matrix to multiply (B) by and the other is the output matrix AB.
CSRMatrix<T>* matMatMult(CSRMatrix<T>& mat_right);

This method transposes a CSR Matrix, but only if it is lower triangular.
void transposeiflower();

Cholesky decomposition decomposes the positive definite CSR matrix into a lower triangular and its conjugate component.
CSRMatrix<T>* CholeskyDecomp();

Example:

    int rows = 3;
    int cols = 3;
    int nnzs = 0;

    auto* A = new Matrix<double>(rows, cols, 3 * pow(rows, 2), 2 * pow(rows, 2));
    auto* sparse_mat = new CSRMatrix<double>(rows, cols, nnzs, true);
    auto* output = new CSRMatrix<double>(rows, cols, nnzs, true);
    sparse_mat->dense2sparse(*A, sparse_mat);
    output = sparse_mat->CholeskyDecomp();

    A->printMatrix();
    output->printMatrix();
    -------------
    Printing matrix

    32 1 1
    1 44 1
    1 1 34
    Printing matrix
    Values: 5.65685 0.176777 6.63089 0.176777 0.146096 5.82644
    row_position: 0 1 3 6
    col_index: 0 0 1 0 1 2

Forward and backward substitutions to solve LUx = b, where U is the upper triangular with units on the diagonal and L is the lower triangular.
void fsubstitution(T* b, T* y);
void bsubstitution(T* y, T* x);

This is a direct solver method using LL<<sup>T</sup>> decomposition with back/forward substitution.
void CholeskySolve(Matrix<T>* A, T* b, T* x);

We hope you enjoyed the read.
