#include <iostream>
#include <math.h>
#include <ctime>
#include "Matrix.h"
#include "Matrix.cpp"
#include "CSRMatrix.h"
#include "CSRMatrix.cpp"

using namespace std;


int main_22()
{

   int  rows = 500;
   int cols = rows;

    // Testing our matrix class
   auto* dense_mat = new Matrix<double>(rows, cols, true);


   int nnzs = 8;
   auto *sparse_mat = new CSRMatrix<double>(rows, cols, nnzs, true);
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

   sparse_mat_2->jacobi_solver_sparse(foo, output, 10, false);




   dense_mat->matVecMult(output, answer_check12);
   RMS_a = dense_mat->RMS_norm_diff(foo, answer_check12);
   cout << "RMS " << RMS_a << endl;


   
   delete[] foo;
   delete[] output;
   delete[] answer_check12;
   delete[] answer_check;
   delete dense_mat;
   delete sparse_mat;
   delete sparse_mat_2;


    system("pause");

    return 0;

}