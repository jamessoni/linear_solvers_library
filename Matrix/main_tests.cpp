#include <iostream>
#include <math.h>
#include <ctime>
#include <vector>
#include <chrono>
#include <memory>
#include <fstream>
#include <string>
#include <sstream> //to split strings
#include <iterator> //to split strings as well
#include "Matrix.h"
#include "Matrix.cpp"
#include "CSRMatrix.h"
#include "CSRMatrix.cpp"

using namespace std;

void test_printValues()
{   
    cout << "printValues():" << "\n";
    int rows = 4;
    int cols = 4;
    auto *dense_mat = new Matrix<double>(rows, cols, true);

    vector<double> vs = {5, 6, 2, 1, 9, 9, 7, 2, 4, 3, 8, 1, 2, 0, 9, 1};
    for (int i = 0; i < rows * cols; i++)
    {
      dense_mat->values[i] = vs[i];
    }
    
    dense_mat->printMatrix();
    delete dense_mat;
}

void test_SPDMatrixcheck()
{
    int rows = 3;
    int cols = 3;

    //creating the SPD matrix
    auto* A1 = new Matrix<double>(rows, cols, 3 * pow(rows, 2), 2 * pow(rows, 2));

    cout << "SPDMatrixcheck(): " <<endl;
    
    A1->printMatrix();
    A1->SPDMatrixcheck(); 
}


void readMatrixFromFile(string name, Matrix<double> *toread)
{
    vector<string> *splitl;
    fstream myfile;
    string line;
    myfile.open(name);
    while (getline(myfile, line)) {
        istringstream iss(line);
        vector<string> splitline((istream_iterator<string>(iss)), istream_iterator<string>());
        for(int i = 0; i<splitline.size();i++)
        {
            toread->values[i] = stod(splitline[i]);
        }
    }
    myfile.close();
}

void test_matMatMult()
{   
    int rows = 100;
    int cols = 100;
    auto *dense_mat1 = new Matrix<double>(rows, cols, true);
    auto *dense_mat2 = new Matrix<double>(rows, cols, true);
    auto *dense_mat3 = new Matrix<double>(rows, cols, true);
    auto *dense_mat4 = new Matrix<double>(rows, cols, true);
    readMatrixFromFile("MMM100-1.txt",dense_mat1);
    readMatrixFromFile("MMM100-2.txt",dense_mat2);
    readMatrixFromFile("MMM100-3.txt",dense_mat3);
    
    clock_t start = clock();
    dense_mat1->matMatMult(*dense_mat2,*dense_mat4);
    
    clock_t end = clock();

    double RMS = dense_mat1->RMS_norm_diff(dense_mat3->values, dense_mat4->values);
    if (RMS > 1.e-2)
    {
        cout << "MatMatMult failed.";
    } else
    {
        cout << "MatMatMult successful.";
    }
    cout << "Time spent to multiply two " << rows <<"x" << cols<< " matrices: " << (double) (end-start) / (double)(CLOCKS_PER_SEC) * 1000.0 << endl;
    
    // dense_mat3->printMatrix();
    // dense_mat4->printMatrix();

    delete dense_mat1;
    delete dense_mat2;
    delete dense_mat3;
    delete dense_mat4;
}


// test_gauss_seidel(Matrix<double>& a, Matrix<double>& b, Matrix<double>& x)
// {
//     void gauss_seidel(a, b, x);
//     a.matVecMult(b);
// }

//    void vecVecsubtract(T* vec_a, T* vec_b, T* output);
//    float RMS_norm_diff(T* vec_a, T* vec_b);

//    // Jacobi solver element-wise
//    void jacobi_solver_element(T* b, T* output, int maxIter, bool initialised);

//    // Functions for Jacobi solver matrix
//    void jacobi_solver_matrix(double* b, double* output, int maxIter, bool initialised);
//    void jacobi_decomposition(Matrix<T>* D, Matrix<T>* N);
	
//    void LUDecomp(Matrix<T>& L, Matrix<T>& U);
//    void SLUDecomp(Matrix<T>& LU);
//    void IPLUDecomp();
   
//    void fsubstitution(Matrix<T>& L, T* y,T* b);
//    void bsubstitution(Matrix<T>& U, T* x, T* y);

//    void LUSolve(double* b, double* output, bool inplace);
//    void conjugate_gradient(T* b, T* x, int maxIter, float tol);



// void test_LUSolve()
// {
//     int rows = 10;
// 	int cols = 10;
// 	auto* A = new Matrix<double>(rows, cols, true);
//     double* b = new double[rows];

// 	for (int i = 0; i < rows; i++) {
// 		b[i] = rand() % 300 + 5;
// 	}

//     readMatrixFromFile("MMM10-1.txt",A);
// 	// set up arrays to store solutions and answer check
// 	double* x = new double[rows];
//     double* answer_check = new double[rows];
//     A->LUSolve(b,x,false);
//     A->matVecMult(x,answer_check);

//     // for (int i = 0; i<10; i++)
//     // {
//     //     cout << b[i] << " " << answer_check[i] << endl;
//     // }

//     // if this is close to zero, the function works 
// 	double RMS = A->RMS_norm_diff(b, answer_check);
// 	// if RMS is larger than e-6, function is not working
// 	if (RMS > 1.e-6) {
// 		cout << " LUSolve failed" << endl;
// 	} else
//     {
//         cout << " LUSolve successful" << endl;
//     }
// 	delete[] b;
// 	delete[] x;
// 	delete[] answer_check;
// 	delete[] A;
// }

void test_jacobi_solver_matrix()
{
    float tol = 1e-6;
    vector<int> vs = {10,100,1000};
    for(int i=0;i<3;i++)
    {
        int rows = vs[i];
	    int cols = vs[i];
        auto *A = new Matrix<double>(rows, cols, rows*rows, cols*cols);

        double* b = new double[rows];

        for (int i = 0; i < rows; i++) {
            b[i] = rand() % 300 + 5;
        }

        // set up arrays to store solutions and answer check
        double* x = new double[rows];
        double* answer_check = new double[rows];
        for(int i = 0; i<rows;i++)
        {
            x[i] = 0.0;
        }
        

        clock_t start = clock();
        A->jacobi_solver_matrix(b,x,1000,true, tol);
    
        clock_t end = clock();

        A->matVecMult(x,answer_check);

        //if this is close to zero, the function works 
        double RMS = A->RMS_norm_diff(b, answer_check);

        // if RMS is larger than e-6, function is not working
        if (RMS > 1.e-2) {
            cout << "Jacobi_solver_matrix method failed. " << "Time spent to solve: " << (double) (end-start) / (double)(CLOCKS_PER_SEC) * 1000.0 << endl;
        } else
        {
            cout << "Jacobi_solver_matrix method successful. "<< "Time spent to solve a "<< rows << "x" << rows <<" matrix: "<< (double) (end-start) / (double)(CLOCKS_PER_SEC) * 1000.0 << endl;
        }
        delete[] b;
        delete[] x;
        delete[] answer_check;
        delete A;
    }
}

void test_jacobi_solver_element()
{
    float tol = 1e-6;
    vector<int> vs = {10,100,1000};
    for(int i=0;i<3;i++)
    {
        int rows = vs[i];
	    int cols = vs[i];
        auto *A = new Matrix<double>(rows, cols, rows*rows, cols*cols);

        double* b = new double[rows];

        for (int i = 0; i < rows; i++) {
            b[i] = rand() % 300 + 5;
        }

        // set up arrays to store solutions and answer check
        double* x = new double[rows];
        double* answer_check = new double[rows];
        for(int i = 0; i<rows;i++)
        {
            x[i] = 0.0;
        }
        

        clock_t start = clock();
        A->jacobi_solver_element(b,x,1000,true, tol);
    
        clock_t end = clock();

        A->matVecMult(x,answer_check);

        //if this is close to zero, the function works 
        double RMS = A->RMS_norm_diff(b, answer_check);

        // if RMS is larger than e-6, function is not working
        if (RMS > 1.e-2) {
            cout << "Jacobi_solver_element method failed. " << "Time spent to solve: " << (double) (end-start) / (double)(CLOCKS_PER_SEC) * 1000.0 << endl;
        } else
        {
            cout << "Jacobi_solver_element method successful. "<< "Time spent to solve a "<< rows << "x" << rows <<" matrix: "<< (double) (end-start) / (double)(CLOCKS_PER_SEC) * 1000.0 << endl;
        }
        delete[] b;
        delete[] x;
        delete[] answer_check;
        delete A;
    }
}

//void test_LUSolve()
//{
//    vector<int> vs = {10};
//    for(int i=0;i<1;i++)
//    {
//        int rows = vs[i];
//	    int cols = vs[i];
//        auto *A = new Matrix<double>(rows, cols, rows*rows, cols*cols);
//        double* b = new double[rows];
//
//        for (int i = 0; i < rows; i++) {
//            b[i] = rand() % 300 + 5;
//        }
//
//        // set up arrays to store solutions and answer check
//        double* x = new double[rows];
//        double* answer_check = new double[rows];
//        for(int i = 0; i<rows;i++)
//        {
//            x[i] = 0.0;
//        }
//        
//
//        clock_t start = clock();
//        A->LUSolve(b,x,false);
//    
//        clock_t end = clock();
//
//        A->matVecMult(x,answer_check);
//
//        //if this is close to zero, the function works 
//        double RMS = A->RMS_norm_diff(b, answer_check);
//
//        // if RMS is larger than e-6, function is not working
//        if (RMS > 1.e-2) {
//            cout << "LUSolve method failed. " << "Time spent to solve: " << (double) (end-start) / (double)(CLOCKS_PER_SEC) * 1000.0 << endl;
//        } else
//        {
//            cout << "LUsolve method successful. "<< "Time spent to solve a "<< rows << "x" << rows <<" matrix: "<< (double) (end-start) / (double)(CLOCKS_PER_SEC) * 1000.0 << endl;
//        }
//        delete b;
//        delete x;
//        delete answer_check;
//        delete A;
//    }
//}


void test_conjugate_gradient()
{
    vector<int> vs = {10,100,1000};
    for(int i=0;i<3;i++)
    {
        int rows = vs[i];
	    int cols = vs[i];
        auto *A = new Matrix<double>(rows, cols, rows*rows, cols*cols);

        double* b = new double[rows];

        for (int i = 0; i < rows; i++) {
            b[i] = rand() % 300 + 5;
        }

        // set up arrays to store solutions and answer check
        double* x = new double[rows];
        double* answer_check = new double[rows];
        for(int i = 0; i<rows;i++)
        {
            x[i] = 0.0;
        }
        
        //add tolerances. For LU use machine precision.

        double tol = 1.e-10;
        clock_t start = clock();
        A->conjugate_gradient(b,x,1000,tol);
    
        clock_t end = clock();

        A->matVecMult(x,answer_check);

        //if this is close to zero, the function works 
        double RMS = A->RMS_norm_diff(b, answer_check);

        // if RMS is larger than e-6, function is not working
        if (RMS > 1.e-2) {
            cout << "Conjugate_gradient method failed. " << "Time spent to solve: " << (double) (end-start) / (double)(CLOCKS_PER_SEC) * 1000.0 << endl;
        } else
        {
            cout << "Conjugate_gradient method successful. "<< "Time spent to solve a "<< rows << "x" << rows <<" matrix: "<< (double) (end-start) / (double)(CLOCKS_PER_SEC) * 1000.0 << endl;
        }
        delete[] b;
        delete[] x;
        delete[] answer_check;
        delete A;
    }
}

void test_gauss_seidel()
{   
    float tol = 1e-6;
    vector<int> vs = {10,100,1000};
    for(int i=0;i<3;i++)
    {
        int rows = vs[i];
        int cols = vs[i];

        //creating the SPD matrix
        auto* A = new Matrix<double>(rows, cols, 3 * pow(rows, 1.5), 2 * pow(rows, 1.5));

        //creating space for solution matrix (x) and matrix (b)

        double* b = new double[rows];
        for (int i = 0; i < rows; i++) {
            b[i] = rand() % 300 + 5;
        }

        // set up arrays to store solutions and answer check
        double* x = new double[rows];
        double* answer_check = new double[rows];
        for (int i = 0; i < rows; i++)
        {
            x[i] = 0.0;
        }

        //using gauss-seidel method
        clock_t start = clock();
        A->gauss_seidel(*A, b, x, tol);
        clock_t end = clock();


        A->matVecMult(x,answer_check);
        double RMS = A->RMS_norm_diff(b, answer_check);

        if (RMS > 1e-5) {
            cout << "Gauss_Seidel method failed. " << "Time spent to solve: " << (double) (end-start) / (double)(CLOCKS_PER_SEC) * 1000.0 << endl;
        } else
        {
            cout << "Gauss_Seidel method successful. "<< "Time spent to solve a "<< rows << "x" << rows <<" matrix: "<< (double) (end-start) / (double)(CLOCKS_PER_SEC) * 1000.0 << endl;
        }

        delete[] b;
        delete[] x;
        delete[] answer_check;
        delete A;
    }
}

//void test_gauss_seidel_sparse()
//{
//    float tol = 1e-6;
//    vector<int> vs = { 10,100,1000 };
//    for (int i = 0; i < 3; i++)
//    {
//        int rows = vs[i];
//        int cols = vs[i];
//
//        //creating dense mat
//        auto* dense_mat = new Matrix<double>(rows, cols, 3 * pow(rows, 1.5), 2 * pow(rows, 1.5));
//
//        int nnzs = 8;
//        auto* sparse_mat = new CSRMatrix<double>(rows, cols, nnzs, true);
//        auto* sparse_mat_2 = new CSRMatrix<double>(rows, cols, nnzs, true);
//        sparse_mat->row_position[nnzs] = nnzs;
//
//        sparse_mat->dense2sparse(*dense_mat, sparse_mat_2);
//
//        //creating space for solution matrix (x) and matrix (b)
//
//        double* b = new double[rows];
//        for (int i = 0; i < rows; i++) {
//            b[i] = rand() % 300 + 5;
//        }
//
//        // set up arrays to store solutions and answer check
//        double* x = new double[rows];
//        double* answer_check = new double[rows];
//        for (int i = 0; i < rows; i++)
//        {
//            x[i] = 0.0;
//        }
//
//        //using gauss-seidel method
//        clock_t start = clock();
//        A->gauss_seidel(*A, b, x, tol);
//        clock_t end = clock();
//
//
//        A->matVecMult(x, answer_check);
//        double RMS = A->RMS_norm_diff(b, answer_check);
//
//        if (RMS > 1e-5) {
//            cout << "Gauss_Seidel method failed. " << "Time spent to solve: " << (double)(end - start) / (double)(CLOCKS_PER_SEC) * 1000.0 << endl;
//        }
//        else
//        {
//            cout << "Gauss_Seidel method successful. " << "Time spent to solve a " << rows << "x" << rows << " matrix: " << (double)(end - start) / (double)(CLOCKS_PER_SEC) * 1000.0 << endl;
//        }
//
//        delete[] b;
//        delete[] x;
//        delete[] answer_check;
//        delete A;
//    }
//
//
//}

//known solution
// x = (-3/65, 1/12, 3/13, 4/13)
void test_sparse_gauss_seidel()
{
    float tol = 1e-6;
    int rows = 4;
    int cols = 4;

    int const nnzs = 7;
    auto* sparse_mat = new CSRMatrix<double>(rows, cols, nnzs, true);

    //sparse_mat->row_position[nnzs] = nnzs;
    int vals[nnzs] = { 10,2,12,6,2,1,9 };
    for (int i = 0; i < nnzs; i++)
    {
        sparse_mat->values[i] = vals[i];
    }

    int col_ind[nnzs] = { 0,2,1,2,3,2,3 };
    for (int i = 0; i < nnzs; i++)
    {
        sparse_mat->col_index[i] = col_ind[i];
    }

    int row_ind[5] = { 0,2,3,5,7 };

    for (int i = 0; i < 5; i++)
    {
        sparse_mat->row_position[i] = row_ind[i];
    }

    // set up arrays to store b, solution and answer check
    double* x = new double[rows];
    double* b = new double[rows];
    double* answer_check = new double[rows];

    //creating space for solution matrix (x) and matrix (b)

    // filling x with 0's as initial guess
    for (int i = 0; i < rows; i++)
    {
        x[i] = 0.0;
    }

    // filling b with i
    for (int i = 0; i < rows; i++)
    {
        b[i] = i;
    }

    //using gauss-seidel method
    clock_t start = clock();
    sparse_mat->gauss_seidel_sparse(*sparse_mat, b, x, tol);
    clock_t end = clock();


    sparse_mat->matVecMult(x, answer_check);
    double RMS = sparse_mat->RMS_norm_diff(b, answer_check);

    if (RMS > 1e-5) {
        cout << "Sparse Gauss Seidel with known input/output method failed. " << "Time spent to solve: " << (double)(end - start) / (double)(CLOCKS_PER_SEC) * 1000.0 << endl;
    }
    else
    {
        cout << "Sparse Gauss Seidel with known input/output method successful. " << "Time spent to solve a " << rows << "x" << rows << " matrix: " << (double)(end - start) / (double)(CLOCKS_PER_SEC) * 1000.0 << endl;
    }

    delete[] b;
    delete[] x;
    delete[] answer_check;
    delete sparse_mat;

}

//known solution
// x = (-3/65, 1/12, 3/13, 4/13)
void test_sparse_jacobi()
{
    float tol = 1e-6;
    int rows = 4;
    int cols = 4;

    int const nnzs = 7;
    auto* sparse_mat = new CSRMatrix<double>(rows, cols, nnzs, true);

    //sparse_mat->row_position[nnzs] = nnzs;
    int vals[nnzs] = { 10,2,12,6,2,1,9 };
    for (int i = 0; i < nnzs; i++)
    {
        sparse_mat->values[i] = vals[i];
    }

    int col_ind[nnzs] = { 0,2,1,2,3,2,3 };
    for (int i = 0; i < nnzs; i++)
    {
        sparse_mat->col_index[i] = col_ind[i];
    }

    int row_ind[5] = { 0,2,3,5,7 };

    for (int i = 0; i < 5; i++)
    {
        sparse_mat->row_position[i] = row_ind[i];
    }

    // set up arrays to store b, solution and answer check
    double* x = new double[rows];
    double* b = new double[rows];
    double* answer_check = new double[rows];

    //creating space for solution matrix (x) and matrix (b)

    // filling x with 0's as initial guess
    for (int i = 0; i < rows; i++)
    {
        x[i] = 0.0;
    }

    // filling b with i
    for (int i = 0; i < rows; i++)
    {
        b[i] = i;
    }

    //using jacobi method
    clock_t start = clock();
    sparse_mat->jacobi_solver_sparse(b, x, 1000, true, tol);
    clock_t end = clock();


    sparse_mat->matVecMult(x, answer_check);
    double RMS = sparse_mat->RMS_norm_diff(b, answer_check);

    if (RMS > 1e-5) {
        cout << "Sparse Jacobi with known input/output method failed. " << "Time spent to solve: " << (double)(end - start) / (double)(CLOCKS_PER_SEC) * 1000.0 << endl;
    }
    else
    {
        cout << "Sparse Jacobi with known input/output method successful. " << "Time spent to solve a " << rows << "x" << rows << " matrix: " << (double)(end - start) / (double)(CLOCKS_PER_SEC) * 1000.0 << endl;
    }

    delete[] b;
    delete[] x;
    delete[] answer_check;
    delete sparse_mat;

}

int main_22_1()
{
    float tol = 1e-6;
    int  rows = 500;
    int cols = rows;

    // Testing our matrix class
    auto* dense_mat = new Matrix<double>(rows, cols, true);


    int nnzs = 8;
    auto* sparse_mat = new CSRMatrix<double>(rows, cols, nnzs, true);
    auto* sparse_mat_2 = new CSRMatrix<double>(rows, cols, nnzs, true);
    sparse_mat->row_position[nnzs] = nnzs;


    sparse_mat->dense2sparse(*dense_mat, sparse_mat_2);

    double* foo = new double[rows];
    for (int i = 0; i < rows; i++) {
        foo[i] = i;
    }
    // set up arrays to store solutions for output to screen
    //double output1[rows];
    double RMS_a = 0;
    double* output = new double[rows];
    double* answer_check = new double[rows];
    double* answer_check12 = new double[rows];

    clock_t start = clock();
    sparse_mat_2->jacobi_solver_sparse(foo, output, 10, false, tol);
    clock_t end = clock();

    dense_mat->matVecMult(output, answer_check);
    double RMS = dense_mat->RMS_norm_diff(foo, answer_check);

    if (RMS > 1e-5) {
        cout << "Sparse Jacobi method failed. " << "Time spent to solve: " << (double)(end - start) / (double)(CLOCKS_PER_SEC) * 1000.0 << endl;
    }
    else
    {
        cout << "Sparse Jacobi method successful. " << "Time spent to solve a " << rows << "x" << rows << " matrix: " << (double)(end - start) / (double)(CLOCKS_PER_SEC) * 1000.0 << endl;
    }

    delete[] foo;
    delete[] output;
    delete[] answer_check12;
    delete[] answer_check;
    delete dense_mat;
    delete sparse_mat;
    delete sparse_mat_2;

    return 0;

}

int Sparse_gauss_seidel_trial()
{
    float tol = 1e-6;
    int  rows = 500;
    int cols = rows;

    // Testing our matrix class
    auto* dense_mat = new Matrix<double>(rows, cols, true);


    int nnzs = 8;
    auto* sparse_mat = new CSRMatrix<double>(rows, cols, nnzs, true);
    auto* sparse_mat_2 = new CSRMatrix<double>(rows, cols, nnzs, true);
    sparse_mat->row_position[nnzs] = nnzs;


    sparse_mat->dense2sparse(*dense_mat, sparse_mat_2);

    double* foo = new double[rows];
    for (int i = 0; i < rows; i++) {
        foo[i] = i;
    }
    // set up arrays to store solutions for output to screen
    //double output1[rows];
    double RMS_a = 0;
    double* output = new double[rows];
    double* answer_check = new double[rows];
    double* answer_check12 = new double[rows];

    clock_t start = clock();
    sparse_mat_2->gauss_seidel_sparse(*sparse_mat, foo, output, tol);
    clock_t end = clock();

    dense_mat->matVecMult(output, answer_check);
    double RMS = dense_mat->RMS_norm_diff(foo, answer_check);

    if (RMS > 1e-5) {
        cout << "Sparse Gauss-seidel method failed. " << "Time spent to solve: " << (double)(end - start) / (double)(CLOCKS_PER_SEC) * 1000.0 << endl;
    }
    else
    {
        cout << "Sparse Gauss-seidel method successful. " << "Time spent to solve a " << rows << "x" << rows << " matrix: " << (double)(end - start) / (double)(CLOCKS_PER_SEC) * 1000.0 << endl;
    }

    delete[] foo;
    delete[] output;
    delete[] answer_check12;
    delete[] answer_check;
    delete dense_mat;
    delete sparse_mat;
    delete sparse_mat_2;

    return 0;

}

// random CSR matrices of increasing size

void test_gauss_seidel_sparse()
{
    float tol = 1.e-6;
    vector<int> vs = { 10,100,1000 };
    for (int i = 0; i < 3; i++)
    {
        int nnzs = 0;
        int rows = vs[i];
        int cols = vs[i];
        auto* A = new Matrix<double>(rows, cols, rows * rows, cols * cols);
        //A->printMatrix();
        
        for (int j = 0; j < A->rows * A->cols; j++)
        {
            int val = A->values[j];
            if (val != 0)
            {
                nnzs += 1;
            }
        }
    
        auto* sparse_mat = new CSRMatrix<double>(rows, cols, nnzs, true);
        auto* sparse_mat_2 = new CSRMatrix<double>(rows, cols, nnzs, true);
        sparse_mat->row_position[nnzs] = nnzs;

        sparse_mat->dense2sparse(*A, sparse_mat_2);

        double* b = new double[rows];

        for (int i = 0; i < rows; i++) {
            b[i] = rand() % 300 + 5;
        }

        // set up arrays to store solutions and answer check
        double* x = new double[rows];
        double* answer_check = new double[rows];
        for (int i = 0; i < rows; i++)
        {
            x[i] = 0.0;
        }

        //add tolerances. For LU use machine precision.

        clock_t start = clock();
        sparse_mat_2->gauss_seidel_sparse(*sparse_mat_2, b, x, tol);

        clock_t end = clock();

        A->matVecMult(x, answer_check);

        //if this is close to zero, the function works 
        double RMS = A->RMS_norm_diff(b, answer_check);

        // if RMS is larger than e-6, function is not working
        if (RMS > 1.e-5) {
            cout << "Sparse Gauss-seidel method failed. " << "Time spent to solve: " << (double)(end - start) / (double)(CLOCKS_PER_SEC) * 1000.0 << endl;
        }
        else
        {
            cout << "Sparse Gauss-seidel method successful. " << "Time spent to solve a " << rows << "x" << rows << " matrix: " << (double)(end - start) / (double)(CLOCKS_PER_SEC) * 1000.0 << endl;
        }
        delete[] b;
        delete[] x;
        delete[] answer_check;
        delete A;
        delete sparse_mat;
        delete sparse_mat_2;
    }

}

void test_jacobi_sparse()
{
    float tol = 1.e-6;
    vector<int> vs = { 10,100,1000 };
    for (int i = 0; i < 3; i++)
    {
        int nnzs = 0;
        int rows = vs[i];
        int cols = vs[i];
        auto* A = new Matrix<double>(rows, cols, rows * rows, cols * cols);
        //A->printMatrix();

        for (int j = 0; j < A->rows * A->cols; j++)
        {
            int val = A->values[j];
            if (val != 0)
            {
                nnzs += 1;
            }
        }

        auto* sparse_mat = new CSRMatrix<double>(rows, cols, nnzs, true);
        auto* sparse_mat_2 = new CSRMatrix<double>(rows, cols, nnzs, true);
        sparse_mat->row_position[nnzs] = nnzs;

        sparse_mat->dense2sparse(*A, sparse_mat_2);

        double* b = new double[rows];

        for (int i = 0; i < rows; i++) {
            b[i] = rand() % 300 + 5;
        }

        // set up arrays to store solutions and answer check
        double* x = new double[rows];
        double* answer_check = new double[rows];
        for (int i = 0; i < rows; i++)
        {
            x[i] = 0.0;
        }

        //add tolerances. For LU use machine precision.

        clock_t start = clock();
        sparse_mat_2->jacobi_solver_sparse(b, x, 100, true, tol);

        clock_t end = clock();

        A->matVecMult(x, answer_check);

        //if this is close to zero, the function works 
        double RMS = A->RMS_norm_diff(b, answer_check);

        // if RMS is larger than e-6, function is not working
        if (RMS > 1.e-5) {
            cout << "Sparse Jacobi method failed. " << "Time spent to solve: " << (double)(end - start) / (double)(CLOCKS_PER_SEC) * 1000.0 << endl;
        }
        else
        {
            cout << "Sparse Jacobi method successful. " << "Time spent to solve a " << rows << "x" << rows << " matrix: " << (double)(end - start) / (double)(CLOCKS_PER_SEC) * 1000.0 << endl;
        }
        delete[] b;
        delete[] x;
        delete[] answer_check;
        delete A;
        delete sparse_mat;
        delete sparse_mat_2;
    }

}

int main()
{
    cout << "Testing:" << "\n\n";
    //test_printValues();
    //test_matMatMult();
    //cout << "\nSPD Dense solver testing: \n";
    //test_jacobi_solver_matrix();
    //test_jacobi_solver_element();
    ////test_LUSolve();
    //test_conjugate_gradient();
    //test_gauss_seidel();

    //sparse solver tests
    cout << "\nSparse solver testing:\n";
    test_gauss_seidel_sparse(); //10x10 100x100 1000x1000
    test_sparse_gauss_seidel(); //known solution
    test_sparse_jacobi(); //known solution
    test_jacobi_sparse(); //10x10 100x100 1000x1000

    //~
    //main_22_1();
    //Sparse_gauss_seidel_trial();


    return 0;
}