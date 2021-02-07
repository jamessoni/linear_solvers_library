#include <iostream>
#include <string>
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
#include <vector>

using namespace std;
using namespace std::chrono;

void test_jacobi_solver_matrix()
{
    float tol = 1e-6;
    vector<int> vs = { 10,100,1000 };
    for (int i = 0; i < 3; i++)
    {
        int rowscols = vs[i];

        auto* A = new Matrix<double>(rowscols, rowscols, 3 * pow(rowscols, 2), 2 * pow(rowscols, 2));

        // set up arrays to store solutions and answer check
        double* x = new double[rowscols];
        double* b = new double[rowscols];
        double* answer_check = new double[rowscols];

        for (int i = 0; i < rowscols; i++) {
            b[i] = rand() % 300 + 5;
        }


        for (int i = 0; i < rowscols; i++)
        {
            x[i] = 0.0;
        }


        auto start = high_resolution_clock::now();
        A->jacobi_solver_matrix(A, b, x, 1000, true, tol);
        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<microseconds>(stop - start);

        A->matVecMult(x, answer_check);

        //if this is close to zero, the function works 
        double RMS = A->RMS_norm_diff(b, answer_check);

        // if RMS is larger than e-6, function is not working
        if (RMS > 1.e-4) {
            cout << "Jacobi solver matrix method failed. " << "Time spent to solve: " << duration.count() << endl;
        }
        else
        {
            cout << "Jacobi solver matrix method successful. " << "Time spent to solve: " << duration.count() << endl;
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
    vector<int> vs = { 10,100,1000 };
    for (int i = 0; i < 3; i++)
    {
        int rows = vs[i];
        int cols = vs[i];
        //auto *A = new Matrix<double>(rows, cols, rows*rows, cols*cols);
        auto* A = new Matrix<double>(rows, cols, 3 * pow(rows, 2), 2 * pow(rows, 2));
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


        auto start = high_resolution_clock::now();
        A->jacobi_solver_element(A, b, x, 1000, true, tol);
        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<microseconds>(stop - start);

        A->matVecMult(x, answer_check);

        //if this is close to zero, the function works 
        double RMS = A->RMS_norm_diff(b, answer_check);

        // if RMS is larger than e-6, function is not working
        if (RMS > 1.e-4) {
            cout << "Jacobi solver element method failed. " << "Time spent to solve: " << duration.count() << endl;
        }
        else
        {
            cout << "Jacobi solver element method successful. " << "Time spent to solve: " << duration.count() << endl;
        }
        delete[] b;
        delete[] x;
        delete[] answer_check;
        delete A;
    }
}

void test_LUSolve()
{
    vector<int> vs = { 10,100,1000 };
    for (int i = 0; i < 3; i++)
    {
        int rows = vs[i];
        int cols = vs[i];
        auto* A = new Matrix<double>(rows, cols, rows * rows, cols * cols);
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

        auto start = high_resolution_clock::now();
        A->LUSolve(A, b, x, false);
        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<microseconds>(stop - start);

        A->matVecMult(x, answer_check);
        //if this is close to zero, the function works 
        double RMS = A->RMS_norm_diff(b, answer_check);

        // if RMS is larger than e-6, function is not working
        if (RMS > 1.e-2) {
            cout << "LU solver method failed. " << "Time spent to solve: " << duration.count() << endl;
        }
        else
        {
            cout << "LU solver method successful. " << "Time spent to solve: " << duration.count() << endl;
        }
        delete[] b;
        delete[] x;
        delete[] answer_check;
        delete A;
    }
}


void test_conjugate_gradient()
{
    double tol = 1.e-6;
    vector<int> vs = { 10,100,1000 };
    for (int i = 0; i < 3; i++)
    {
        int rows = vs[i];
        int cols = vs[i];
        //auto *A = new Matrix<double>(rows, cols, rows*rows, cols*cols);
        auto* A = new Matrix<double>(rows, cols, 3 * pow(rows, 2), 2 * pow(rows, 2));
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

        auto start = high_resolution_clock::now();
        A->conjugate_gradient(A, b, x, 1000, tol);
        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<microseconds>(stop - start);

        A->matVecMult(x, answer_check);

        //if this is close to zero, the function works 
        double RMS = A->RMS_norm_diff(b, answer_check);

        // if RMS is larger than e-6, function is not working
        if (RMS > 1.e-4) {
            cout << "Conjugate gradient solver method failed. " << "Time spent to solve: " << duration.count() << endl;
        }
        else
        {
            cout << "Conjugate gradient solver method successful. " << "Time spent to solve: " << duration.count() << endl;
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
    vector<int> vs = { 10,100,1000 };
    for (int i = 0; i < 3; i++)
    {
        int rows = vs[i];
        int cols = vs[i];

        //creating the SPD matrix
        auto* A = new Matrix<double>(rows, cols, 3 * pow(rows, 1.5), 2 * pow(rows, 1.5));

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
        auto start = high_resolution_clock::now();
        A->gauss_seidel(A, b, x, 1000, tol);
        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<microseconds>(stop - start);


        A->matVecMult(x, answer_check);
        double RMS = A->RMS_norm_diff(b, answer_check);

        if (RMS > 1e-4) {
            cout << "Gauss-seidel solver method failed. " << "Time spent to solve: " << duration.count() << endl;
        }
        else
        {
            cout << "Gauss-seidel solver method successful. " << "Time spent to solve: " << duration.count() << endl;
        }

        delete[] b;
        delete[] x;
        delete[] answer_check;
        delete A;
    }
}


//known solution
void test_gauss_seidel_sparse_known()
{
    int max_Iter = 500;
    float tol = 1e-6;
    int rows = 4;
    int cols = 4;

    int const nnzs = 8;
    auto* sparse_mat = new CSRMatrix<double>(rows, cols, nnzs, true);

    //sparse_mat->row_position[nnzs] = nnzs;
    //int vals[nnzs] = { 10,2,12,6,2,1,9 };
    int vals[nnzs] = { 10,2,12,2,6,1,1,9 };
    for (int i = 0; i < nnzs; i++)
    {
        sparse_mat->values[i] = vals[i];
    }

    int col_ind[nnzs] = { 0,2,1,0,2,3,2,3 };
    for (int i = 0; i < nnzs; i++)
    {
        sparse_mat->col_index[i] = col_ind[i];
    }

    int row_ind[5] = { 0,2,3,6,8 };

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
    auto start = high_resolution_clock::now();
    sparse_mat->gauss_seidel_sparse(sparse_mat, b, x, max_Iter, tol);
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);


    sparse_mat->matVecMult(x, answer_check);
    double RMS = sparse_mat->RMS_norm_diff(b, answer_check);

    if (RMS > 1e-4) {
        cout << "Gauss-seidel sparse solver method failed for known input-output. " << "Time spent to solve: " << duration.count() << endl;
    }
    else
    {
        cout << "Gauss-seidel sparse solver method successful for known input-output. " << "Time spent to solve: " << duration.count() << endl;
    }

    delete[] b;
    delete[] x;
    delete[] answer_check;
    delete sparse_mat;

}

//known solution
void test_jacobi_sparse_known()
{
    float tol = 1e-6;
    int rows = 4;
    int cols = 4;

    int const nnzs = 8;
    auto* sparse_mat = new CSRMatrix<double>(rows, cols, nnzs, true);

    //sparse_mat->row_position[nnzs] = nnzs;
    int vals[nnzs] = { 10,2,12,2,6,1,1,9 };

    for (int i = 0; i < nnzs; i++)
    {
        sparse_mat->values[i] = vals[i];
    }

    int col_ind[nnzs] = { 0,2,1,0,2,3,2,3 };
    for (int i = 0; i < nnzs; i++)
    {
        sparse_mat->col_index[i] = col_ind[i];
    }

    /* int row_ind[5] = { 0,2,3,5,7 };*/
    int row_ind[5] = { 0,2,3,6,8 };
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
    auto start = high_resolution_clock::now();
    sparse_mat->jacobi_solver_sparse(sparse_mat, b, x, 1000, true, tol);
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);


    sparse_mat->matVecMult(x, answer_check);
    double RMS = sparse_mat->RMS_norm_diff(b, answer_check);

    if (RMS > 1e-4) {
        cout << "Jacobi sparse solver method failed for known input-output. " << "Time spent to solve: " << duration.count() << endl;
    }
    else
    {
        cout << "Jacobi sparse solver method successful for known input-output. " << "Time spent to solve: " << duration.count() << endl;
    }

    delete[] b;
    delete[] x;
    delete[] answer_check;
    delete sparse_mat;

}

//known solution
void test_conjugate_gradient_known()
{
    float tol = 1e-6;
    int rows = 4;
    int cols = 4;

    int const nnzs = 8;
    auto* sparse_mat = new CSRMatrix<double>(rows, cols, nnzs, true);

    //sparse_mat->row_position[nnzs] = nnzs;
    int vals[nnzs] = { 10,2,12,2,6,1,1,9 };

    for (int i = 0; i < nnzs; i++)
    {
        sparse_mat->values[i] = vals[i];
    }

    int col_ind[nnzs] = { 0,2,1,0,2,3,2,3 };
    for (int i = 0; i < nnzs; i++)
    {
        sparse_mat->col_index[i] = col_ind[i];
    }

    /* int row_ind[5] = { 0,2,3,5,7 };*/
    int row_ind[5] = { 0,2,3,6,8 };
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

    //using CG method
    auto start = high_resolution_clock::now();
    sparse_mat->conjugate_gradient(sparse_mat, b, x, 1000, tol);
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);

    sparse_mat->matVecMult(x, answer_check);
    double RMS = sparse_mat->RMS_norm_diff(b, answer_check);
    
    if (RMS > 1e-4) {
        cout << "Conjugate gradient sparse solver method failed for known input-output. " << "Time spent to solve: " << duration.count() << endl;
    }
    else
    {
        cout << "Conjugate gradient sparse solver method failed for known input-output. " << "Time spent to solve: " << duration.count() << endl;
    }

    delete[] b;
    delete[] x;
    delete[] answer_check;
    delete sparse_mat;
}

//known solution
void test_cholesky_known()
{
    float tol = 1e-6;
    int rows = 4;
    int cols = 4;

    int const nnzs = 8;
    auto* sparse_mat = new CSRMatrix<double>(rows, cols, nnzs, true);

    //sparse_mat->row_position[nnzs] = nnzs;
    int vals[nnzs] = { 10,2,12,2,6,1,1,9 };

    for (int i = 0; i < nnzs; i++)
    {
        sparse_mat->values[i] = vals[i];
    }

    int col_ind[nnzs] = { 0,2,1,0,2,3,2,3 };
    for (int i = 0; i < nnzs; i++)
    {
        sparse_mat->col_index[i] = col_ind[i];
    }

    /* int row_ind[5] = { 0,2,3,5,7 };*/
    int row_ind[5] = { 0,2,3,6,8 };
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

    //using cholesky method
    auto start = high_resolution_clock::now();
    sparse_mat->CholeskySolve(sparse_mat, b, x);
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);

    sparse_mat->matVecMult(x, answer_check);
    double RMS = sparse_mat->RMS_norm_diff(b, answer_check);
    
    if (RMS > 1e-4) {
        cout << "Cholesky sparse solver method failed for known input-output. " << "Time spent to solve: " << duration.count() << endl;
    }
    else
    {
        cout << "Cholesky sparse solver method failed for known input-output. " << "Time spent to solve: " << duration.count() << endl;
    }

    delete[] b;
    delete[] x;
    delete[] answer_check;
    delete sparse_mat;
}

// random CSR matrices of increasing size
void test_gauss_seidel_sparse()
{
    int max_iter = 500;
    float tol = 1.e-6;
    vector<int> vs = { 10,100,1000 };
    //vector<int> vs = { 10,100 };
    for (int i = 0; i < 3; i++)
    {
        int nnzs = 0;
        int rows = vs[i];
        int cols = vs[i];
        auto* A = new Matrix<double>(rows, cols, 3 * pow(rows, 2), 2 * pow(rows, 2));

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

        auto start = high_resolution_clock::now();
        sparse_mat_2->gauss_seidel_sparse(sparse_mat_2, b, x, max_iter, tol);
        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<microseconds>(stop - start);

        A->matVecMult(x, answer_check);

        //if this is close to zero, the function works 
        double RMS = A->RMS_norm_diff(b, answer_check);

        // if RMS is larger than e-6, function is not working
        if (RMS > 1.e-4) {
            cout << "Gauss-seidel sparse solver method failed. " << "Time spent to solve: " << duration.count() << endl;
        }
        else
        {
            cout << "Gauss-seidel sparse solver method successful. " << "Time spent to solve: " << duration.count() << endl;
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
        //auto* A = new Matrix<double>(rows, cols, rows * rows, cols * cols);
        auto* A = new Matrix<double>(rows, cols, 3 * pow(rows, 2), 2 * pow(rows, 2));
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

        auto start = high_resolution_clock::now();
        sparse_mat_2->jacobi_solver_sparse(sparse_mat_2, b, x, 100, true, tol);
        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<microseconds>(stop - start);

        A->matVecMult(x, answer_check);

        //if this is close to zero, the function works 
        double RMS = A->RMS_norm_diff(b, answer_check);

        // if RMS is larger than e-6, function is not working
        if (RMS > 1.e-4) {
            cout << "Jacobi sparse solver method failed. " << "Time spent to solve: " << duration.count() << endl;
        }
        else
        {
            cout << "Jacobi sparse solver method successful. " << "Time spent to solve: " << duration.count() << endl;
        }
        delete[] b;
        delete[] x;
        delete[] answer_check;
        delete A;
        delete sparse_mat;
        delete sparse_mat_2;
    }

}

void test_choleskyDecomp()
{
    vector<int> vs = { 10,100,1000 };
    for (int i = 0; i < vs.size(); i++)
    {
        int rowscols = vs[i];

        auto* A = new Matrix<double>(rowscols, rowscols, 3 * pow(rowscols, 1.5), 2 * pow(rowscols, 1.5));

        double* b = new double[rowscols];
        double* x = new double[rowscols];

        // filling x with 0's as initial guess
        // filling b with random ints
        for (int i = 0; i < rowscols; i++)
        {
            b[i] = rand() % 10;
            x[i] = 0;
        }

        auto* L = new Matrix<double>(rowscols, rowscols, true);

        auto start = high_resolution_clock::now();
        A->CholeskySolve(A, b, x);
        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<microseconds>(stop - start);

        double* answer_check = new double[rowscols];
        A->matVecMult(x, answer_check);

        double RMS = A->RMS_norm_diff(b, answer_check);

        if (RMS > 1.e-2) {
            cout << "Cholesky solver method failed. " << "Time spent to solve: " << duration.count() << endl;
        }
        else
        {
            cout << "Cholesky solver method successful. " << "Time spent to solve: " << duration.count() << endl;
        }


        delete A;
        delete L;
        delete[] x;
        delete[] b;
        delete[] answer_check;
    }
}

void test_sparse_CholeskySolve()
{
    vector<int> vs = { 10,100 };
    for (int i = 0; i < vs.size(); i++)
    {
        int rowscols = vs[i];

        auto* A = new Matrix<double>(rowscols, rowscols, 3 * pow(rowscols, 1.5), 2 * pow(rowscols, 1.5));


        double* b = new double[rowscols];
        double* x = new double[rowscols];

        // filling x with 0's as initial guess
        // filling b with random ints
        for (int i = 0; i < rowscols; i++)
        {
            b[i] = rand() % 10;
            x[i] = 0;
        }

        int nnzs = pow(rowscols, 2);
        auto* sparse_mat = new CSRMatrix<double>(rowscols, rowscols, nnzs, true);

        sparse_mat->dense2sparse(*A, sparse_mat);

        auto start = high_resolution_clock::now();
        sparse_mat->CholeskySolve(sparse_mat, b, x);
        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<microseconds>(stop - start);

        double* answer_check = new double[rowscols];
        A->matVecMult(x, answer_check);

        double RMS = A->RMS_norm_diff(b, answer_check);

        if (RMS > 1.e-2) {
            cout << "Cholesky sparse solver method failed. " << "Time spent to solve: " << duration.count() << endl;
        }
        else
        {
            cout << "Cholesky sparse solver method successful. " << "Time spent to solve: " << duration.count() << endl;
        }

        delete[] b;
        delete[] x;
        delete[] answer_check;
        delete A;
        delete sparse_mat;
    }
}

void test_sparse_conjugate_gradient()
{
    vector<int> vs = { 10,100,1000 };
    for (int i = 0; i < vs.size(); i++)
    {
        double tol = 1.e-6;
        int rowscols = vs[i];

        auto* A = new Matrix<double>(rowscols, rowscols, 3 * pow(rowscols, 1.5), 2 * pow(rowscols, 1.5));


        double* b = new double[rowscols];
        double* x = new double[rowscols];

        // filling x with 0's as initial guess
        // filling b with random ints
        for (int i = 0; i < rowscols; i++)
        {
            b[i] = rand() % 10;
            x[i] = 0;
        }

        int nnzs = pow(rowscols, 2);
        auto* sparse_mat = new CSRMatrix<double>(rowscols, rowscols, nnzs, true);

        sparse_mat->dense2sparse(*A, sparse_mat);

        auto start = high_resolution_clock::now();
        sparse_mat->conjugate_gradient(sparse_mat, b, x, 1000, tol);
        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<microseconds>(stop - start);

        double* answer_check = new double[rowscols];
        A->matVecMult(x, answer_check);

        double RMS = A->RMS_norm_diff(b, answer_check);

        if (RMS > 1.e-2) {
            cout << "Conjugate gradient sparse solver method failed. " << "Time spent to solve: " << duration.count() << endl;
        }
        else
        {
            cout << "Conjugate gradient sparse solver method successful. " << "Time spent to solve: " << duration.count() << endl;
        }

        delete[] b;
        delete[] x;
        delete[] answer_check;
        delete A;
        delete sparse_mat;
    }
}

void gauss_seidel_example()
{

    int rows = 0;
    int cols = 0;
    float tol = 1e-6;

    cout << "Input the number of Matrix rows: ";
    cin >> rows;
    //ensuring a square matrix
    cols = rows;

    //creating the SPD matrix
    auto* A = new Matrix<double>(rows, cols, 3 * pow(rows, 2), 2 * pow(rows, 2));

    A->printMatrix();

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

    cout << "initialised input x: ";
    for (int i = 0; i < rows; i++)
    {
        cout << x[i] << " ";
    }
    cout << endl;

    cout << "randomly initialised input b: ";
    for (int i = 0; i < rows; i++)
    {
        cout << b[i] << " ";
    }
    cout << endl;

    A->gauss_seidel(A, b, x, 1000, tol);

    cout << "solution x: ";
    for (int i = 0; i < rows; i++)
    {
        cout << x[i] << " ";
    }
    cout << endl << endl;

    delete[] b;
    delete[] x;
    delete[] answer_check;
    delete A;

}

void jacobi_matrix_example()
{
    int rows = 0;
    int cols = 0;
    float tol = 1e-6;

    cout << "Input the number of Matrix rows: ";
    cin >> rows;
    //ensuring a square matrix
    cols = rows;

    //creating the SPD matrix
    auto* A = new Matrix<double>(rows, cols, 3 * pow(rows, 2), 2 * pow(rows, 2));

    A->printMatrix();

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

    cout << "initialised input x: ";
    for (int i = 0; i < rows; i++)
    {
        cout << x[i] << " ";
    }
    cout << endl;

    cout << "randomly initialised input b: ";
    for (int i = 0; i < rows; i++)
    {
        cout << b[i] << " ";
    }
    cout << endl;

    A->jacobi_solver_matrix(A, b, x, 100, true, tol);

    cout << "solution x: ";
    for (int i = 0; i < rows; i++)
    {
        cout << x[i] << " ";
    }
    cout << endl << endl;

    delete[] b;
    delete[] x;
    delete[] answer_check;
    delete A;

}

void jacobi_element_example()
{
    int rows = 0;
    int cols = 0;
    float tol = 1e-6;

    cout << "Input the number of Matrix rows: ";
    cin >> rows;
    //ensuring a square matrix
    cols = rows;

    //creating the SPD matrix
    auto* A = new Matrix<double>(rows, cols, 3 * pow(rows, 2), 2 * pow(rows, 2));

    A->printMatrix();

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

    cout << "initialised input x: ";
    for (int i = 0; i < rows; i++)
    {
        cout << x[i] << " ";
    }
    cout << endl;

    cout << "randomly initialised input b: ";
    for (int i = 0; i < rows; i++)
    {
        cout << b[i] << " ";
    }
    cout << endl;

    A->jacobi_solver_element(A, b, x, 100, true, tol);

    cout << "solution x: ";
    for (int i = 0; i < rows; i++)
    {
        cout << x[i] << " ";
    }
    cout << endl << endl;

    delete[] b;
    delete[] x;
    delete[] answer_check;
    delete A;
}

void LU_example()
{
    int rows = 0;
    int cols = 0;
    float tol = 1e-6;

    cout << "Input the number of Matrix rows: ";
    cin >> rows;
    //ensuring a square matrix
    cols = rows;

    //creating the SPD matrix
    auto* A = new Matrix<double>(rows, cols, 3 * pow(rows, 2), 2 * pow(rows, 2));

    A->printMatrix();

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

    cout << "initialised input x: ";
    for (int i = 0; i < rows; i++)
    {
        cout << x[i] << " ";
    }
    cout << endl;

    cout << "randomly initialised input b: ";
    for (int i = 0; i < rows; i++)
    {
        cout << b[i] << " ";
    }
    cout << endl;

    A->LUSolve(A, b, x, false);

    cout << "solution x: ";
    for (int i = 0; i < rows; i++)
    {
        cout << x[i] << " ";
    }
    cout << endl << endl;

    delete[] b;
    delete[] x;
    delete[] answer_check;
    delete A;
}

void conjugate_gradient_example()
{
    int rows = 0;
    int cols = 0;
    float tol = 1e-6;

    cout << "Input the number of Matrix rows: ";
    cin >> rows;
    //ensuring a square matrix
    cols = rows;

    //creating the SPD matrix
    auto* A = new Matrix<double>(rows, cols, 3 * pow(rows, 2), 2 * pow(rows, 2));

    A->printMatrix();

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

    cout << "initialised input x: ";
    for (int i = 0; i < rows; i++)
    {
        cout << x[i] << " ";
    }
    cout << endl;

    cout << "randomly initialised input b: ";
    for (int i = 0; i < rows; i++)
    {
        cout << b[i] << " ";
    }
    cout << endl;

    A->conjugate_gradient(A, b, x, 100, tol);

    cout << "solution x: ";
    for (int i = 0; i < rows; i++)
    {
        cout << x[i] << " ";
    }
    cout << endl << endl;

    delete[] b;
    delete[] x;
    delete[] answer_check;
    delete A;
}

void cholesky_example()
{
    int rows = 0;
    int cols = 0;
    float tol = 1e-6;

    cout << "Input the number of Matrix rows: ";
    cin >> rows;
    //ensuring a square matrix
    cols = rows;

    //creating the SPD matrix
    auto* A = new Matrix<double>(rows, cols, 3 * pow(rows, 2), 2 * pow(rows, 2));

    A->printMatrix();

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

    cout << "initialised input x: ";
    for (int i = 0; i < rows; i++)
    {
        cout << x[i] << " ";
    }
    cout << endl;

    cout << "randomly initialised input b: ";
    for (int i = 0; i < rows; i++)
    {
        cout << b[i] << " ";
    }
    cout << endl;

    A->CholeskySolve(A, b, x);

    cout << "solution x: ";
    for (int i = 0; i < rows; i++)
    {
        cout << x[i] << " ";
    }
    cout << endl << endl;

    delete[] b;
    delete[] x;
    delete[] answer_check;
    delete A;
}

void jacobi_sparse_example()
{
    int nnzs = 0;
    int rows = 0;
    int cols = 0;
    float tol = 1e-6;

    cout << "Input the number of Matrix rows: ";
    cin >> rows;
    //ensuring a square matrix
    cols = rows;

    //creating the SPD matrix
    auto* A = new Matrix<double>(rows, cols, 3 * pow(rows, 2), 2 * pow(rows, 2));

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

    sparse_mat_2->jacobi_solver_sparse(sparse_mat_2, b, x, 100, true, tol);

    cout << "initialised input x: ";
    for (int i = 0; i < rows; i++)
    {
        cout << x[i] << " ";
    }
    cout << endl;

    cout << "randomly initialised input b: ";
    for (int i = 0; i < rows; i++)
    {
        cout << b[i] << " ";
    }
    cout << endl;


    cout << "solution x: ";
    for (int i = 0; i < rows; i++)
    {
        cout << x[i] << " ";
    }
    cout << endl << endl;

    delete[] b;
    delete[] x;
    delete[] answer_check;
    delete A;
    delete sparse_mat;
    delete sparse_mat_2;
}

void Gauss_seidel_sparse_example()
{
    int max_iter = 500;
    int nnzs = 0;
    int rows = 0;
    int cols = 0;
    float tol = 1e-6;

    cout << "Input the number of Matrix rows: ";
    cin >> rows;
    //ensuring a square matrix
    cols = rows;

    //creating the SPD matrix
    auto* A = new Matrix<double>(rows, cols, 3 * pow(rows, 2), 2 * pow(rows, 2));

    A->printMatrix();

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

    sparse_mat_2->gauss_seidel_sparse(sparse_mat_2, b, x, max_iter, tol);

    cout << "initialised input x: ";
    for (int i = 0; i < rows; i++)
    {
        cout << x[i] << " ";
    }
    cout << endl;

    cout << "randomly initialised input b: ";
    for (int i = 0; i < rows; i++)
    {
        cout << b[i] << " ";
    }
    cout << endl;


    cout << "solution x: ";
    for (int i = 0; i < rows; i++)
    {
        cout << x[i] << " ";
    }
    cout << endl << endl;

    delete[] b;
    delete[] x;
    delete[] answer_check;
    delete A;
    delete sparse_mat;
    delete sparse_mat_2;

}

void conjugant_gradient_sparse_example()
{
    int nnzs = 0;
    int rows = 0;
    int cols = 0;
    float tol = 1e-6;

    cout << "Input the number of Matrix rows: ";
    cin >> rows;
    //ensuring a square matrix
    cols = rows;

    //creating the SPD matrix
    auto* A = new Matrix<double>(rows, cols, 3 * pow(rows, 2), 2 * pow(rows, 2));

    A->printMatrix();

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

    sparse_mat_2->conjugate_gradient(sparse_mat_2, b, x, 1000, tol);

    cout << "initialised input x: ";
    for (int i = 0; i < rows; i++)
    {
        cout << x[i] << " ";
    }
    cout << endl;

    cout << "randomly initialised input b: ";
    for (int i = 0; i < rows; i++)
    {
        cout << b[i] << " ";
    }
    cout << endl;


    cout << "solution x: ";
    for (int i = 0; i < rows; i++)
    {
        cout << x[i] << " ";
    }
    cout << endl << endl;

    delete[] b;
    delete[] x;
    delete[] answer_check;
    delete A;
    delete sparse_mat;
    delete sparse_mat_2;
}

void cholesky_sparse_example()
{
    int nnzs = 0;
    int rows = 0;
    int cols = 0;
    float tol = 1e-6;

    cout << "Input the number of Matrix rows: ";
    cin >> rows;
    //ensuring a square matrix
    cols = rows;

    //creating the SPD matrix
    auto* A = new Matrix<double>(rows, cols, 3 * pow(rows, 2), 2 * pow(rows, 2));

    A->printMatrix();

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

    sparse_mat_2->CholeskySolve(sparse_mat_2, b, x);

    cout << "initialised input x: ";
    for (int i = 0; i < rows; i++)
    {
        cout << x[i] << " ";
    }
    cout << endl;

    cout << "randomly initialised input b: ";
    for (int i = 0; i < rows; i++)
    {
        cout << b[i] << " ";
    }
    cout << endl;


    cout << "solution x: ";
    for (int i = 0; i < rows; i++)
    {
        cout << x[i] << " ";
    }
    cout << endl << endl;

    delete[] b;
    delete[] x;
    delete[] answer_check;
    delete A;
    delete sparse_mat;
    delete sparse_mat_2;

}

int main_test_solvers()
{
    cout << "Testing in microseconds:" << "\n";
    cout << "\nSPD Dense solver testing: \n";
    test_jacobi_solver_matrix();
    test_jacobi_solver_element();
    test_LUSolve();
    test_conjugate_gradient();
    test_gauss_seidel();
    test_choleskyDecomp();

    //sparse solver tests
    cout << "\nSparse solver testing:\n";
    test_gauss_seidel_sparse(); //10x10 100x100 1000x1000
    test_jacobi_sparse(); //10x10 100x100 1000x1000
    test_sparse_conjugate_gradient(); //10x10 100x100 1000x1000
    test_sparse_CholeskySolve(); //10x10 100x100

    cout << "\nKnown input-output testing: \n";
    test_gauss_seidel_sparse_known(); //known solution
    test_jacobi_sparse_known(); //known solution
    //test_conjugate_gradient_known(); //known solution
    //test_cholesky_known(); //known solution
    cout << endl << endl;
    return 0;
}
void UI()
{
    char input_str;
    cout << "Team Void's Linear solver library\n\n";
    cout << "Do you want to solve try a dense or sparse solver?\n"
        "Type 'D' for dense or 'S' for sparse or 'T' for solver convergence comparison tests or 'E' to exit:\n"
        "choice: ";
    cin >> input_str;

    if (input_str == 'd' || input_str == 'D')
    {
        int solver_choice;
        cout << "Dense solvers: ";
        cout << "\n1. Jacobi matrix solver \n2. Jacobi element solver\n3. Gauss-seidel Solver\n4. LU Solver\n5. Conjugate Gradient Solver\n6. Cholesky Solver\n";
        cout << "Type the number allocated to the solver you wish to use: ";
        cin >> solver_choice;

        if (solver_choice == 1)
        {
            jacobi_matrix_example();
            UI();
        }
        else if (solver_choice == 2)
        {
            jacobi_element_example();
            UI();
        }
        else if (solver_choice == 3)
        {
            gauss_seidel_example();
            UI();
        }
        else if (solver_choice == 4)
        {
            LU_example();
            UI();
        }
        else if (solver_choice == 5)
        {
            conjugate_gradient_example();
            UI();
        }
        else if (solver_choice == 6)
        {
            cholesky_example();
            UI();
        }
    }
    else if (input_str == 's' || input_str == 'S')
    {
        int solver_choice;
        cout << "Sparse solvers: ";
        cout << "\n1. Jacobi sparse solver \n2. Gauss-seidel sparse solver\n3. Conjugant gradient sparse solver\n4. Cholesky sparse solver\n";
        cout << "Type the number allocated to the solver you wish to use: ";
        cin >> solver_choice;

        if (solver_choice == 1)
        {
            jacobi_sparse_example();
            UI();
        }
        if (solver_choice == 2)
        {
            Gauss_seidel_sparse_example();
            UI();
        }
        if (solver_choice == 3)
        {
            conjugant_gradient_sparse_example();
            UI();
        }
        if (solver_choice == 4)
        {
            cholesky_sparse_example();
            UI();
        }
    }
    else if (input_str == 't' || input_str == 'T')
    {
        main_test_solvers();
        UI();
    }
    else if (input_str == 'e' || input_str == 'E')
    {
        cout << "\nHave a great day!\n";
        system("pause");
        exit(3);
    }
}

int main()
{
    UI();

    return 0;
}