#include <iostream>
#include <math.h>
#include <ctime>
#include "Matrix.h"
#include "Matrix.cpp"

using namespace std;

int main()
{

   int rows = 100;
   int cols = 100;

   // Testing our matrix class
   auto *dense_mat = new Matrix<double>(rows, cols, true);


   // Now let's test printing our matrix with our other function
   //dense_mat->printMatrix();      
   //ddense_mat->printMatrix();



   // create b vector and fill it with random values
   double foo[100];
   for (int i = 0; i < 100; i++) {
	   foo[i] = rand() % 300 + 5;
   }
   // set up arrays to store solutions for output to screen
   double output1[100];
   double* output = output1;{}
   double* answer_check = new double[100];
   double* answer_check12 = new double[100];

 


   clock_t start = clock();
   dense_mat->jacobi_solver_matrix(foo, output, 100000);
   clock_t end = clock();
   cout << "Time spent on matrix version: " << (double)(end - start) / (double)(CLOCKS_PER_SEC) * 1000.0 << endl;

   clock_t start2 = clock();
   dense_mat->jacobi_solver_element(foo, output, 100000);
   clock_t end2 = clock();
   cout << "Time spent element-wise version: " << (double)(end2 - start2) / (double)(CLOCKS_PER_SEC) * 1000.0 << endl;



   // uncomment below check solution, if output to screen is close to zero than last solver is working correct


   /*
   dense_mat->matVecMult(output, answer_check12);

   dense_mat->vecVecsubtract(foo, answer_check12, answer_check);

   
   for (int i = 0; i < dense_mat->cols; i++) {
	   cout << answer_check[i] << endl;
   }
   
   */
   

   delete dense_mat;
   system("pause");

}