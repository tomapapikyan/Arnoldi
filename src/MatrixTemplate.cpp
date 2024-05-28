#include "MatrixTemplate.h"
#include <iostream>

Matvec::Matvec(int n, int m){
    matrix = new double[n*m];
    columns = m;
    rows = n;
}

Matvec::Matvec(Matvec const& other):
    columns(other.columns), rows(other.rows), matrix(new double[rows*columns]){
    for (int i=0; i<rows*columns; i++)
        matrix[i] = other.matrix[i];

}

double*
Matvec::get_matrix(){
    return matrix;
}

int
Matvec::get_row_count(){
    return rows;
}

int
Matvec::get_column_count(){
    return columns;
}

void
Matvec::mat_vec(double* vec, double* prod){
omp_set_num_threads(6);
#pragma omp parallel for collapse(2) reduction(+:prod[:rows])
for (int i=0; i<rows; i++)
    for (int j=0; j<columns; j++)
        prod[i] += matrix[j*rows + i] * vec[j];
}

Matvec::~Matvec(){
    delete[] matrix;
}

void Matvec::fill_matrix(){
for (int i=0; i<rows; i++)
    for (int j=0; j<=i; j++){
        matrix[j*rows+i] = (j+1)%4000;
        matrix[i*rows+j] = (j+1)%4000;
    }
}
