#include <iostream>
#include <fstream>
#include "Arnoldi_Iterations.h"
#include "CppBLAS.h"

int main(){
    std::ifstream input("../../test/input.txt");
    double *A, *wr, *wi;
    int n, m;

    if (!input.is_open())
        return 1;

    input>>n;
    A = new double[n*n];
    for (int i=0; i<n; i++)
        for (int j=0; j<n; j++)
            input>>A[j*n+i];
    std::cout<<"Number of eigenvalues needed:"<<std::endl;
    std::cin>>m;
    wr = new double[2*m];//wr, wi - real and complex parts of eigenvalues
    wi = new double[2*m];
    Compute_Eigenvalues(n, m, A, wr, wi);
    delete[] A;
    delete[] wr;
    delete[] wi;

}
