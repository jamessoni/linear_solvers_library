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

float tol = 1.e-2;
vector<int> vs = {10,100};


void test_jacobi_solver_matrix()
{
    for(int i=0;i<vs.size();i++)
    {
        int rowscols = vs[i];
        
        auto* A = new Matrix<double>(rowscols, rowscols, 3 * pow(rowscols,   1.10), 2 * pow(rowscols,   1.10));
        // set up arrays to store solutions and answer check
        double* x = new double[rowscols];
        double* b = new double[rowscols];
        double* answer_check = new double[rowscols];
        //filling b with rands and x with 0
        for (int i = 0; i < rowscols; i++) {
            b[i] = rand() % 300 + 5;
            x[i] = 0.0;
        }
        

        clock_t start = clock();
        A->jacobi_solver_matrix(A,b,x,1000,true, tol);
    
        clock_t end = clock();

        A->matVecMult(x,answer_check);

        //if this is close to zero, the function works 
        double RMS = A->RMS_norm_diff(b, answer_check);
         
        // if RMS is larger than e-6, function is not working
        if (RMS > tol) {
            cout << "Jacobi_solver_matrix method failed. " << "Time spent to solve: " << (double) (end-start) / (double)(CLOCKS_PER_SEC) * 1000.0 << endl;
        } else
        {
            cout << "Jacobi_solver_matrix method successful. "<< "Time spent to solve a "<< rowscols << "x" << rowscols <<" matrix: "<< (double) (end-start) / (double)(CLOCKS_PER_SEC) * 1000.0 << endl;
        }
        delete[] b;
        delete[] x;
        delete[] answer_check;
        delete A;
    }
}

void test_jacobi_solver_element()
{
    for(int i=0;i<vs.size();i++)
    {
        int rows = vs[i];
	    int cols = vs[i];
        //auto *A = new Matrix<double>(rows, cols, rows*rows, cols*cols);
        auto* A = new Matrix<double>(rows, cols, 3 * pow(rows,   1.10), 2 * pow(rows,   1.10));
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
        A->jacobi_solver_element(A,b,x,1000,true, tol);
    
        clock_t end = clock();

        A->matVecMult(x,answer_check);
        
        //if this is close to zero, the function works 
        double RMS = A->RMS_norm_diff(b, answer_check);
         
        // if RMS is larger than e-6, function is not working
        if (RMS > tol) {
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

void test_LUSolve()
{
    for(int i=0;i<vs.size()-1;i++)
    {
        int rows = vs[i];
	    int cols = vs[i];
        auto *A = new Matrix<double>(rows, cols, 3* pow(rows,  1.10), 2* pow(rows,  1.10));
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
        A->LUSolve(A,b,x,false);
        clock_t end = clock();

        A->matVecMult(x,answer_check);
        //if this is close to zero, the function works 
        double RMS = A->RMS_norm_diff(b, answer_check);
         
        // if RMS is larger than e-6, function is not working
        if (RMS > tol) {
            cout << "LUSolve method failed. " << "Time spent to solve: " << (double) (end-start) / (double)(CLOCKS_PER_SEC) * 1000.0 << endl;
        } else
        {
            cout << "LUsolve method successful. "<< "Time spent to solve a "<< rows << "x" << rows <<" matrix: "<< (double) (end-start) / (double)(CLOCKS_PER_SEC) * 1000.0 << endl;
        }
        delete b;
        delete x;
        delete answer_check;
        delete A;
    }
}


void test_conjugate_gradient()
{
    for(int i=0;i<vs.size();i++)
    {
        int rows = vs[i];
	    int cols = vs[i];
        //auto *A = new Matrix<double>(rows, cols, rows*rows, cols*cols);
        auto* A = new Matrix<double>(rows, cols, 3*pow(rows,  1.10), 2*pow(rows,  1.10));
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

        clock_t start = clock();
        A->conjugate_gradient(A,b,x,1000,tol);
    
        clock_t end = clock();

        A->matVecMult(x,answer_check);

        //if this is close to zero, the function works 
        double RMS = A->RMS_norm_diff(b, answer_check);
         
        // if RMS is larger than e-6, function is not working
        if (RMS > tol) {
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
    for(int i=0;i<vs.size();i++)
    {
        int rows = vs[i];
        int cols = vs[i];

        //creating the SPD matrix
        auto* A = new Matrix<double>(rows, cols, 3 * pow(rows,   1.10), 2 * pow(rows,   1.10));

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
        A->gauss_seidel(A, b, x, 1000, tol);
        clock_t end = clock();


        A->matVecMult(x,answer_check);
        double RMS = A->RMS_norm_diff(b, answer_check);
         
        if (RMS > tol) {
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


//known solution
void test_sparse_gauss_seidel()
{
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

    int col_ind[nnzs] = { 0,2,1,0,2,3,2,3};
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
    clock_t start = clock();
    sparse_mat->gauss_seidel_sparse(sparse_mat, b, x, 50000, tol);
    clock_t end = clock();


    sparse_mat->matVecMult(x, answer_check);
    double RMS = sparse_mat->RMS_norm_diff(b, answer_check);

    if (RMS > tol) {
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
void test_sparse_jacobi()
{
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
    clock_t start = clock();
    sparse_mat->jacobi_solver_sparse(sparse_mat, b, x, 1000, true, tol);
    clock_t end = clock();


    sparse_mat->matVecMult(x, answer_check);
    double RMS = sparse_mat->RMS_norm_diff(b, answer_check);
     
    if (RMS > tol) {
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

// random CSR matrices of increasing size
void test_gauss_seidel_sparse()
{
    for(int i=0;i<vs.size();i++)
    {
        int nnzs = 0;
        int rows = vs[i];
        int cols = vs[i];
        //auto* A = new Matrix<double>(rows, cols, rows * rows, cols * cols);
        auto* A = new Matrix<double>(rows, cols, 3 * pow(rows,   1.10), 2 * pow(rows,   1.10));
        //A->SPDMatrixcheck();
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

        //add tolerances. For LU use machine precision.

        clock_t start = clock();
        sparse_mat_2->gauss_seidel_sparse(sparse_mat_2, b, x, 50000, tol);

        clock_t end = clock();

        A->matVecMult(x, answer_check);

        //if this is close to zero, the function works 
        double RMS = A->RMS_norm_diff(b, answer_check);
         
        // if RMS is larger than e-6, function is not working
        if (RMS > tol) {
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
    for(int i=0;i<vs.size();i++)
    {
        int nnzs = 0;
        int rows = vs[i];
        int cols = vs[i];
        //auto* A = new Matrix<double>(rows, cols, rows * rows, cols * cols);
        auto* A = new Matrix<double>(rows, cols, 3 * pow(rows,   1.10), 2 * pow(rows,   1.10));
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

        //add tolerances. For LU use machine precision.

        clock_t start = clock();
        sparse_mat_2->jacobi_solver_sparse(sparse_mat_2, b, x, 100, true, tol);
        clock_t end = clock();

        A->matVecMult(x, answer_check);

        //if this is close to zero, the function works 
        double RMS = A->RMS_norm_diff(b, answer_check);
         
        // if RMS is larger than e-6, function is not working
        if (RMS > tol) {
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

void test_choleskyDecomp()
{   
    for(int i=0;i<vs.size();i++)
    {
        int rowscols = vs[i];
        
        auto* A = new Matrix<double>(rowscols, rowscols, 3 * pow(rowscols,   1.10), 2 * pow(rowscols,   1.10));

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

        clock_t start = clock();
        A->CholeskySolve(A,b,x);
        clock_t end = clock();

        double* answer_check = new double[rowscols];
        A->matVecMult(x,answer_check);

        double RMS = A->RMS_norm_diff(b, answer_check);

        if (RMS > tol) {
            cout << "CholeskySolve method failed. " << "Time spent to solve: " << (double) (end-start) / (double)(CLOCKS_PER_SEC) * 1000.0 << endl;
        } else
        {
            cout << "CholeskySolve method successful. "<< "Time spent to solve a "<< rowscols << "x" << rowscols <<" matrix: "<< (double) (end-start) / (double)(CLOCKS_PER_SEC) * 1000.0 << endl;
        }

         
        delete A;
        delete L;
        delete x;
        delete b;
        delete answer_check;
    }
}

void test_sparse_CholeskySolve()
{
    for(int i=0;i<vs.size();i++)
    {
        int rowscols = vs[i];
        
        auto* A = new Matrix<double>(rowscols, rowscols, 3 * pow(rowscols,   1.10), 2 * pow(rowscols,   1.10));


        double* b = new double[rowscols];
        double* x = new double[rowscols];

        // filling x with 0's as initial guess
        // filling b with random ints
        for (int i = 0; i < rowscols; i++)
        {
            b[i] = rand() % 10;
            x[i] = 0;
        }

        int nnzs = pow(rowscols,2);
        auto* sparse_mat = new CSRMatrix<double>(rowscols, rowscols, nnzs, true);

        sparse_mat->dense2sparse(*A, sparse_mat);

        clock_t start = clock();
        sparse_mat->CholeskySolve(A,b,x);
        clock_t end = clock();
        double* answer_check = new double[rowscols];
        A->matVecMult(x,answer_check);

        double RMS = A->RMS_norm_diff(b, answer_check);

        if (RMS > tol) {
            cout << "sparse_CholeskySolve method failed. " << "Time spent to solve: " << (double) (end-start) / (double)(CLOCKS_PER_SEC) * 1000.0 << endl;
        } else
        {
            cout << "sparse_CholeskySolve method successful. "<< "Time spent to solve a "<< rowscols << "x" << rowscols <<" matrix: "<< (double) (end-start) / (double)(CLOCKS_PER_SEC) * 1000.0 << endl;
        }
         
        delete A;
        delete sparse_mat;
    }
}

void test_sparse_conjugate_gradient()
{
    for(int i=0;i<vs.size();i++)
    {   
        int rowscols = vs[i];
        
        auto* A = new Matrix<double>(rowscols, rowscols, 3 * pow(rowscols,   1.10), 2 * pow(rowscols,   1.10));


        double* b = new double[rowscols];
        double* x = new double[rowscols];

        // filling x with 0's as initial guess
        // filling b with random ints
        for (int i = 0; i < rowscols; i++)
        {
            b[i] = rand() % 10;
            x[i] = 0;
        }

        int nnzs = pow(rowscols,2);
        auto* sparse_mat = new CSRMatrix<double>(rowscols, rowscols, nnzs, true);

        sparse_mat->dense2sparse(*A, sparse_mat);

        clock_t start = clock();
        sparse_mat->conjugate_gradient(sparse_mat,b,x,1000,tol);
        clock_t end = clock();

        double* answer_check = new double[rowscols];
        A->matVecMult(x,answer_check);

        double RMS = A->RMS_norm_diff(b, answer_check);

          
        if (RMS > tol) {
            cout << "sparse_ConjugateGradient method failed. " << "Time spent to solve: " << (double) (end-start) / (double)(CLOCKS_PER_SEC) * 1000.0 << endl;
        } else
        {
            cout << "sparse_ConjugateGradient method successful. "<< "Time spent to solve a "<< rowscols << "x" << rowscols <<" matrix: "<< (double) (end-start) / (double)(CLOCKS_PER_SEC) * 1000.0 << endl;
        }

        delete A;
        delete sparse_mat;
    }
}

int main()
{
    cout << "Testing:" << "\n\n";
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
    test_sparse_conjugate_gradient();
    test_sparse_CholeskySolve();

    // test_sparse_gauss_seidel(); //known solution
    // test_sparse_jacobi(); //known solution
    //test_sparse_matMatMult();
    return 0;
}