#include <omp.h>

class Matvec{

private:
    int rows;
    int columns;
    double* matrix;

public:
    Matvec(int, int);

    Matvec(Matvec const&);

    ~Matvec();

    double* get_matrix();

    int get_row_count();

    int get_column_count();

    void mat_vec(double*, double*);

    void fill_matrix();
};
