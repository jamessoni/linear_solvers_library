#include <iostream>
#include <math.h>
#include <ctime>
#include "Matrix.h"
#include "Matrix.cpp"

using namespace std;

int main()
{

   int rows = 2;
   int cols = 2;

   // Testing our matrix class
   auto *dense_mat = new Matrix<double>(rows, cols, true);
   // Now we need to go and fill our matrix
   // Let's just fill it with random values for now
   for (int i = 0; i < rows * cols; i++)
   {
      dense_mat->values[i] = i;
   }
   // Now let's test printing our matrix
   dense_mat->printValues();
   // Now let's test printing our matrix with our other function
   dense_mat->printMatrix();   

   // Don't forget to call delete (ie the destructor) once we're done, otherwise 
   // we will leak memory from the new called inside our matrix
   delete dense_mat;

   // ~~~~~~~~~~~~~~~~~~~~~~~
   // ~~~~~~~~~~~~~~~~~~~~~~~
   // What if we allocated the memory ourselves, then decided we want to
   // use it wrapped in our new matrix class
   // The only real difference here is thinking about who "owns" our double array
   // ~~~~~~~~~~~~~~~~~~~~~~~
   // ~~~~~~~~~~~~~~~~~~~~~~~

   auto *matrix_data = new double[rows * cols];
   dense_mat = new Matrix<double>(rows, cols, matrix_data);

   // Now let's fill it, but using the variable matrix_data (which should 
   // be pointing to the same place as dense_mat->values)
   for (int i = 0; i < rows * cols; i++)
   {
      matrix_data[i] = i;
   }
   // Now let's test printing our matrix
   dense_mat->printValues();
   // Now let's test printing our matrix with our other function
   dense_mat->printMatrix();      

   // ~~~~~~~~~~~~~~~~~~~~~~~
   // ~~~~~~~~~~~~~~~~~~~~~~~
   // Now let's test our matmatmult
   // ~~~~~~~~~~~~~~~~~~~~~~~
   // ~~~~~~~~~~~~~~~~~~~~~~~   
   auto *left_mat = new Matrix<double>(rows, cols, matrix_data);
   // Let's make it allocate the space for the output
   auto *output_mat = new Matrix<double>(rows, cols, true);

   // Let's call our matrix multiply routine
   // Let's time it too
   clock_t start = clock();
   dense_mat->matMatMult(*left_mat, *output_mat);
   clock_t end = clock();
   cout << "Time spent: " << (double) (end-start) / (double)(CLOCKS_PER_SEC) * 1000.0 << endl;   

   // Let's check the output
   output_mat->printMatrix();

   //testing the SPD matrix constructor

   //create UI menu method later
   cout << "Team Void's Matrix Library" << endl;
   cout << "Input the number of Matrix rows: ";
   cin >> rows;
   //ensuring a square matrix
   cols = rows;

   //creating a SPD matrix
   auto* spd_mat = new Matrix<double>(rows, cols, 600, 400, 5, 1);

   // Now let's test printing our matrix
   spd_mat->printValues();
   spd_mat->printMatrix();


   // We are now explicitly required to delete our memory
   // As well as deleted our Matrix object
   delete spd_mat;
   delete[] matrix_data;
   delete dense_mat;
   delete left_mat;
   delete output_mat;

	system("pause");

    return 0;
}