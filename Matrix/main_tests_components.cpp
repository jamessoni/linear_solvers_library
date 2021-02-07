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