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
    if (RMS > 1.e-6)
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


void test_jacobi_solver_matrix()
{
    float tol = 1e-6;
    vector<int> vs = {10,100,1000};
    for(int i=0;i<3;i++)
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
        if (RMS > 1.e-4) {
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
        //auto *A = new Matrix<double>(rows, cols, rows*rows, cols*cols);
        auto* A = new Matrix<double>(rows, cols, 3 * pow(rows, 2), 2 * pow(rows, 2));
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
        if (RMS > 1.e-4) {
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
        A->LUSolve(b,x,false);
        clock_t end = clock();

        A->matVecMult(x,answer_check);
        //if this is close to zero, the function works 
        double RMS = A->RMS_norm_diff(b, answer_check);

        // if RMS is larger than e-6, function is not working
        if (RMS > 1.e-2) {
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
    double tol = 1.e-6;
    vector<int> vs = {10,100,1000};
    for(int i=0;i<3;i++)
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
        for(int i = 0; i<rows;i++)
        {
            x[i] = 0.0;
        }
        
        //add tolerances. For LU use machine precision.

        clock_t start = clock();
        A->conjugate_gradient(b,x,1000,tol);
    
        clock_t end = clock();

        A->matVecMult(x,answer_check);

        //if this is close to zero, the function works 
        double RMS = A->RMS_norm_diff(b, answer_check);

        // if RMS is larger than e-6, function is not working
        if (RMS > 1.e-4) {
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

        if (RMS > 1e-4) {
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
    sparse_mat->gauss_seidel_sparse(*sparse_mat, b, x, tol);
    clock_t end = clock();


    sparse_mat->matVecMult(x, answer_check);
    double RMS = sparse_mat->RMS_norm_diff(b, answer_check);

    if (RMS > 1e-4) {
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
    clock_t start = clock();
    sparse_mat->jacobi_solver_sparse(sparse_mat, b, x, 1000, true, tol);
    clock_t end = clock();


    sparse_mat->matVecMult(x, answer_check);
    double RMS = sparse_mat->RMS_norm_diff(b, answer_check);

    if (RMS > 1e-4) {
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
    float tol = 1.e-6;
    vector<int> vs = {10,100,1000 };
    //vector<int> vs = { 10,100 };
    for (int i = 0; i < 3; i++)
    {
        int nnzs = 0;
        int rows = vs[i];
        int cols = vs[i];
        //auto* A = new Matrix<double>(rows, cols, rows * rows, cols * cols);
        auto* A = new Matrix<double>(rows, cols, 3 * pow(rows, 2), 2 * pow(rows, 2));
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
        sparse_mat_2->gauss_seidel_sparse(*sparse_mat_2, b, x, tol);

        clock_t end = clock();

        A->matVecMult(x, answer_check);

        //if this is close to zero, the function works 
        double RMS = A->RMS_norm_diff(b, answer_check);

        // if RMS is larger than e-6, function is not working
        if (RMS > 1.e-4) {
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

        //add tolerances. For LU use machine precision.

        clock_t start = clock();
        sparse_mat_2->jacobi_solver_sparse(sparse_mat_2, b, x, 100, true, tol);
        clock_t end = clock();

        A->matVecMult(x, answer_check);

        //if this is close to zero, the function works 
        double RMS = A->RMS_norm_diff(b, answer_check);

        // if RMS is larger than e-6, function is not working
        if (RMS > 1.e-4) {
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
    vector<int> vs = {10,100,1000};
    for(int i=0;i<3;i++)
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

        clock_t start = clock();
        A->CholeskySolve(b,x);
        clock_t end = clock();

        double* answer_check = new double[rowscols];
        A->matVecMult(x,answer_check);

        double RMS = A->RMS_norm_diff(b, answer_check);

        if (RMS > 1.e-2) {
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

test_sparse_matMatMult()
{
    int rows = 5;
    int cols = 5;
    auto* dense_mat = new Matrix<double>(rows, cols, true);
    auto* dense_mat2 = new Matrix<double>(rows, cols, true);
    auto* dense_mat3 = new Matrix<double>(rows, cols, true);

    vector<double> vs = {0,1,2,3,0,0,0,0,0,0,0,0,7,0,0,8,0,0,0,0,0,0,1,1,1};

    for(int i=0;i<rows*cols;i++)
    {
        dense_mat->values[i] = vs[i];
        dense_mat2->values[i] = vs[i];
    }

    dense_mat->matMatMult(*dense_mat2,*dense_mat3);

    dense_mat->printMatrix();
    dense_mat3->printMatrix();

    int nnzs = 8;
    auto* sparse_mat = new CSRMatrix<double>(rows, cols, nnzs, true);
    auto* sparse_mat2 = new CSRMatrix<double>(rows, cols, nnzs, true);


    sparse_mat->dense2sparse(*dense_mat, sparse_mat);
    sparse_mat->dense2sparse(*dense_mat, sparse_mat2);

    sparse_mat->printMatrix();

    auto* sparse = sparse_mat->matMatMult(*sparse_mat2);

    sparse->printMatrix();

    delete dense_mat;
    delete dense_mat2;
    delete dense_mat3;
    delete sparse_mat;
    delete sparse_mat2;
    delete sparse;
}

test_sparse_Cholesky()
{
    int rows = 5;
    int cols = 5;
    auto* dense_mat = new Matrix<double>(rows, cols, true);
    auto* dense_mat2 = new Matrix<double>(rows, cols, true);

    vector<double> vs = {10, 1, 2, 3, 0,1, 8, 0, 0, 0, 2, 0, 7, 0, 1, 3, 0, 0, 10, 1, 0, 0, 1, 1, 1};

    for(int i=0;i<rows*cols;i++)
    {
        dense_mat->values[i] = vs[i];
    }

    dense_mat->CholeskyDecomp(dense_mat2);

    dense_mat->printMatrix();
    dense_mat2->printMatrix();

    int nnzs = 8;
    auto* sparse_mat = new CSRMatrix<double>(rows, cols, nnzs, true);

    sparse_mat->dense2sparse(*dense_mat, sparse_mat);


    double x[5] = {0,0,0,0,0};
    double b[5] = {1,2,1,2,1};
    double answer_check[5] = {0,0,0,0,0};
    cout << answer_check[0];
    clock_t start = clock();
    sparse_mat->CholeskySolve(b,x);
    clock_t end = clock();
    //cout << "Time taken:" << (double) (end-start) / (double)(CLOCKS_PER_SEC) * 1000.0 << endl;

    sparse_mat->matVecMult(x,answer_check);

    cout << answer_check[0];
    for(int i = 0; i<5; i++)
    {   
        cout << "Heya";
        cout << answer_check[i] << "-"<< b[i] << " ";
    }
    delete dense_mat;
    delete dense_mat2;
    delete sparse_mat;
}

int main()
{
    cout << "Testing:" << "\n\n";
    // test_printValues();
    // test_matMatMult();
    cout << "\nSPD Dense solver testing: \n";
    // test_jacobi_solver_matrix();
    // test_jacobi_solver_element();
    // test_LUSolve();
    // test_conjugate_gradient();
    // test_gauss_seidel();
    // test_choleskyDecomp();

    //sparse solver tests
    cout << "\nSparse solver testing:\n";
    // test_gauss_seidel_sparse(); //10x10 100x100 1000x1000
    // test_sparse_gauss_seidel(); //known solution
    // test_sparse_jacobi(); //known solution
    // test_jacobi_sparse(); //10x10 100x100 1000x1000

    //test_sparse_matMatMult();
    test_sparse_Cholesky();
    return 0;
}