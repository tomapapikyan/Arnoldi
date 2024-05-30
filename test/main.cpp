#include <iostream>
#include <fstream>
#include "Arnoldi_Iterations.hpp"
#include "CppBLAS.h"

int main(){
    int n, m;
    double *wr, *wi, *matrix;
    double itime, ftime, exec_time;

    std::cout<<"Matrix size:\n";
    std::cin>>n;

    matrix = new double[n*n];

    Matvec A(n, n, matrix);
    A.fill_matrix();

    std::cout<<"Number of eigenvalues needed:"<<std::endl;
    std::cin>>m;
    wr = new double[2*m];//wr, wi - real and complex parts of eigenvalues
    wi = new double[2*m];

    itime = omp_get_wtime();
    Compute_Eigenvalues(m, A, wr, wi);
    ftime = omp_get_wtime();

    exec_time = ftime - itime;
    std::cout<<"Time: "<<exec_time<<std::endl;

    delete[] wr;
    delete[] wi;

}
