#include <omp.h>
#include <cmath>
#include <iostream>

#define NUM_THREADS 6

void dgemv(int m, int N, double* A, double* q, double *v){
omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for collapse(2) reduction(+:v[:m])
for (int i=0; i<m; i++)
    for (int j=0; j<N; j++)
        v[i] += A[j*m + i]*q[j];
}


void dgemvT(int m, int N, double* A, double* q, double *v){
omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for collapse(2) reduction(+:v[:m])
for (int i=0; i<m; i++)
    for (int j=0; j<N; j++)
        v[i] += A[i*N + j]*q[j];
}

double ddot(int N, double* q, double* v){
    double res = 0;
    omp_set_num_threads(NUM_THREADS);
    #pragma omp parallel for reduction(+:res)
    for (int i=0; i<N; i++)
        res += q[i]*v[i];
    return res;
}

void daxpy(int N, double a, double* v, double* q){
    omp_set_num_threads(NUM_THREADS);
    #pragma omp parallel for reduction(+:q[:N])
    for (int i=0; i<N; i++)
        q[i] += a*v[i];
}

double dnrm2(int N, double* v){
    double h = 0;
    omp_set_num_threads(NUM_THREADS);
    #pragma omp parallel for reduction(+:h)
    for (int i=0; i<N; i++)
        h += v[i]*v[i];
    return sqrt(h);
}

void dscal(int N, double a, double* v){
    omp_set_num_threads(NUM_THREADS);
    #pragma omp parallel for reduction(*:v[:N])
    for (int i=0; i<N; i++)
        v[i] *= a;
}
