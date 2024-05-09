#include <iostream>
#include "CppBLAS.h"
#include "MatVec.h"
#include <vector>
#include <cstdlib>
#include <omp.h>

#define NUM_THREADS 4

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
        dgemv(N, A, Krylov_bases.back(), v);//v = Aq

        #pragma omp parallel for
        for (int k=0; k< Krylov_bases.size(); k++)//Gram-Schmidt
            Hm[i*m + k] = ddot(N, Krylov_bases[k], v);//h[i,k]=(q[k],v)

        omp_set_num_threads(NUM_THREADS);
        #pragma omp parallel
        {
            double* v_temp = new double[N];
            clear_m(v_temp, N);
            int id = omp_get_thread_num();
            int nthrds = omp_get_num_threads();
            for (int k=id; k<Krylov_bases.size(); k+=nthrds)
                daxpy(N, Hm[i*m + k], Krylov_bases[k], v_temp);//v_temp = v_temp + h*q[j]
            #pragma omp critical
                daxpy(N, -1, v_temp, v);//v = v - v_temp
            delete[] v_temp;
        }

        if (i==(m-1)){
            delete[] v;
            continue;
        }

        double h = dnrm2(N, v);//h=||v||
        Hm[i*m + Krylov_bases.size()] = h;
        dscal(N, 1/h, v);//v = v/h
        Krylov_bases.push_back(v);
    }

    for (double* v: Krylov_bases)
        delete[] v;
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
    delete[] Hm;
}
