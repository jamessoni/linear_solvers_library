#Team Void ACSE-5 Group Assignment

###Implementing Linear Solvers

Within this Advanced Programming Assignment the project aimed to implement an algorithm to solve the linear system *A**x* = *b*, where A is a positive definite matrix, with x and b being vectors.
This linear solver library is written in the programming language, C++, maximising its advantages of Object Oriented Design, memory allocation and optimisation.
The assignment further evaluates, tests and aims to optimise the linear solvers created fulfilling convergence to solutions for both dense and sparse input matrices. Present within this linear solvers library, exists dense and sparse libraries - each with their respective .cpp and .h files. 

####Getting Started / Pre-requisites
#####Mac/Linux:
IDE or C++ compiler of choice. Recommendations listed below:
*default g++ compiler
*Xcode
#####Windows users:
IDE or C++ compiler of choice. Recommendations listed below:
*Microsoft Visual Studio Community
*Visual Studio Code

#####Installation
Clone or download the [Github repository](https://github.com/acse-2020/group-project-team-void) to your local machine

#####Using
For use of the program will depend on your OS and IDE/compiler of choice:
Visual Studio Community IDE users:
*Open a Visual Studio Project, add all the cloned files to the source files list
*Build project and run.

Terminal:
*Enter directory with cloned repository files within the Matrix/ folder
'''
$> g++ [file_name].cpp -fpermissive
$> .\a
'''

####Declaring input Matrices:
#####Dense Matrices:
The values of the dense matrices will depend on whether you as the user wish to input a specified matrix.
One example shown below:
'''
    int rows = 4;
    int cols = 4;

    //creating matrix
    auto *dense_mat = new Matrix<double>(rows, cols, true);

    //allocating specific matrix entries
    vector<double> vs = {5, 6, 2, 1, 9, 9, 7, 2, 4, 3, 8, 1, 2, 0, 9, 1};

    //filling Matrix dense_mat with values specified above
    for (int i = 0; i < rows * cols; i++)
    {
      dense_mat->values[i] = vs[i];
    }

    [...]

    delete dense_mat; 
'''
The example above declares a dense matrix which has further been filled with the values present in the vector (vs). 
However, in providing ease of use of the system - there is functuality to create a random SPD matrix if necessary.
Shown below:
///
