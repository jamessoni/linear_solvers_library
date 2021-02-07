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

void test_printMatrix()
{   
    cout << "Testing printValues():" << "\n";
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
    cout << "\n";
}

void test_SPDMatrixcheck()
{
    int rows = 3;
    int cols = 3;

    //creating the SPD matrix
    auto* A1 = new Matrix<double>(rows, cols, 3 * pow(rows, 2), 2 * pow(rows, 2));

    cout << "Testing SPDMatrixcheck(): " <<endl;
    
    A1->printMatrix();
    A1->SPDMatrixcheck(); 
}


void readMatrixFromFile(string name, Matrix<double> *toread)
{
    vector<string> *splitl;
    fstream myfile;
    string line;
    myfile.open("premade_matrices/"+name);
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
    cout << "Testing matMatMult():" << "\n";
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
    cout << "Time spent to multiply two " << rows <<"x" << cols<< " matrices: " << (double) (end-start) / (double)(CLOCKS_PER_SEC) * 1000.0 << "\n\n";

    delete dense_mat1;
    delete dense_mat2;
    delete dense_mat3;
    delete dense_mat4;
}

void test_transpose()
{
    cout << "Testing transpose():" << "\n";
    int rows = 4;
    int cols = 4;
    auto *dense_mat = new Matrix<double>(rows, cols, true);

    vector<double> vs = {5, 6, 2, 1, 9, 9, 7, 2, 4, 3, 8, 1, 2, 0, 9, 1};
    for (int i = 0; i < rows * cols; i++)
    {
      dense_mat->values[i] = vs[i];
    }
    
    cout << "Input Matrix:" << endl;
    dense_mat->printMatrix();
    dense_mat->transpose();
    cout << "Transposed Matrix:" << endl;
    dense_mat->printMatrix();
    delete dense_mat;
    cout << "\n";
}

void test_matVecMult()
{
    cout << "Testing transpose():" << "\n";
    cout << "\n";
}

void test_sparse_matMatMult()
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

void test_sparse_Cholesky()
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
    cout << "Time taken:" << (double) (end-start) / (double)(CLOCKS_PER_SEC) * 1000.0 << endl;

    sparse_mat->matVecMult(x,answer_check);

    cout << answer_check[0];
    for(int i = 0; i<5; i++)
    {
        cout << answer_check[i] << "-"<< b[i] << " ";
    }
    delete dense_mat;
    delete dense_mat2;
    delete sparse_mat;
}

int main_test_components()
{
    cout << "Testing components:" << "\n\n";
    // test_printMatrix();
    // test_matMatMult();
    // test_transpose();
    // test_matVecMult();

    // test_sparse_matMatMult();
    // test_sparse_Cholesky();
    return 0;
}