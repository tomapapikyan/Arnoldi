#include <iostream>
#include "CppBLAS.h"
#include <vector>
#include <cstdlib>
#include <omp.h>

int compare(const void* a, const void* b){//for qsort
    const double* x = (double*) a;
    const double* y = (double*) b;

    if(*x>*y)
        return -1;
    else if (*x<*y)
        return 1;
    return 0;
}

void clear_m(double* arr, int n){
    for (int i=0; i<n; i++)
        arr[i] = 0;
}

void print_matrix(double* a, int n, int m){
for (int i=0; i<n; i++){
    for (int j=0; j<m; j++)
        std::cout<<a[j*n+i]<<" ";
    std::cout<<std::endl;
    }
}

void Arnoldi_Iteration(double* Hm, double* A, int N, int m){
    double *v = new double[N];
    clear_m(v, N);
    v[0] = 1;

    std::vector<double*> Krylov_bases;
    Krylov_bases.push_back(v);

    for (int i=0; i<m; i++){
        double *v = new double[N];
        clear_m(v, N);
        dgemv_cpp(N, A, Krylov_bases.back(), v);//v = Aq
        #pragma omp parallel for
        for (int k=0; k< Krylov_bases.size(); k++){//Gram-Schmidt
            double* q_j = Krylov_bases[k];
            double h = ddot_cpp(N, q_j, v);//h = (q[j], v)
            Hm[i*m + k] = h;
        }
        for (int k=0; k<Krylov_bases.size(); k++)
            daxpy_cpp(N, -Hm[i*m + k], Krylov_bases[k], v);//v = v - h*q[j]
        if (i==m)
            continue;
        double h = dnrm2_cpp(N, v);//h=||v||
        Hm[i*m + Krylov_bases.size()] = h;//4, 5
        dscal_cpp(N, 1/h, v);//v = v/h
        Krylov_bases.push_back(v);
    }
}

void Compute_Eigenvalues(int N, int k, double* A, double* wr, double* wi){
    int m = 2*k;
    if (m>N)
        m=N;
    double* Hm;
    Hm = new double[m*m];
    clear_m(Hm, m*m);
    Arnoldi_Iteration(Hm, A, N, m);
    dhseqr_cpp(m, Hm, wr, wi); //writes eigenvalues to wr and wi
    std::qsort(wr, m, sizeof(double), compare);
    std::cout<<"First "<<k<<" eigenvalues:"<<std::endl;
    print_matrix(wr, 1, k);
}
