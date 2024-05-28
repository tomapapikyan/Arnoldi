#include <iostream>
#include "CppBLAS.h"
#include "VectorFunctions.h"
#include "Arnoldi_Iterations.h"
#include <vector>
#include <cstdlib>
#include <cmath>

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
    omp_set_num_threads(6);
    #pragma omp parallel
    {
        int id = omp_get_thread_num();
        int nthrds = omp_get_num_threads();
        int k = n/nthrds;
        int begin = id*k;
        int end = begin+k;
        if (id==nthrds-1)
            end = n;
        for (int i=begin; i<end; i++)
            arr[i] = 0;
    }
}

void print_matrix(double* a, int n, int m){
for (int i=0; i<n; i++){
    for (int j=0; j<m; j++)
        std::cout<<a[j*n+i]<<" ";
    std::cout<<std::endl;
    }
}

bool is_null(double* v, int N){
    bool res = true;
    #pragma omp parallel
    {
        int id = omp_get_thread_num();
        int nthrds = omp_get_num_threads();
        int k = N/nthrds;
        int begin = id*k;
        int end = begin+k;
        if (id == nthrds-1)
            end = N;

        for (int i=begin; i<end; i++)
            if(v[i] != 0)
                res = false;
    }
    return res;
}

void Arnoldi_Iteration(double* Hm, Matvec A, int N, int m){

    double* Krylov_bases = new double[m*N];
    clear_m(Krylov_bases, m*N);

    for (int i=0; i<N; i++)
        Krylov_bases[i] = 1/sqrt(N);//q[0]

    double *v = new double[N];
    double *q_temp = new double[N];

    for (int i=0; i<m; i++){
        clear_m(v, N);
        clear_m(q_temp, N);
        A.mat_vec(Krylov_bases + i*N, v); //v = Aq[i]

        dgemvT(i+1, N, Krylov_bases, v, Hm+i*m);//H[i] = Q*q[i]
        dgemv(N, i+1, Krylov_bases, Hm+i*m, q_temp);//q = QH[i]
        daxpy(N, -1, q_temp, v);//v = v - q

        if (i==(m-1))
            continue;

        if (is_null(v, N))
            break;

        double h = dnrm2(N, v);//h=||v||
        Hm[i*m + i+1] = h;
        dscal(N, 1/h, v);//v = v/h
        for (int j=0; j<N; j++)
            Krylov_bases[(i+1)*N + j] = v[j];
    }
    delete[] v;
    delete[] q_temp;
    delete[] Krylov_bases;
}

void Compute_Eigenvalues(int k, Matvec &A, double* wr, double* wi){
    int N = A.get_column_count();
    int m = 2*k;
    if (m>N)
        m=N;
    double* Hm = new double[m*m];
    clear_m(Hm, m*m);
    Arnoldi_Iteration(Hm, A, N, m);
    dhseqr_cpp(m, Hm, wr, wi); //writes eigenvalues to wr and wi
    std::qsort(wr, m, sizeof(double), compare);
//    std::cout<<"First "<<k<<" eigenvalues:"<<std::endl;
//    print_matrix(wr, 1, k);
    delete[] Hm;
}
