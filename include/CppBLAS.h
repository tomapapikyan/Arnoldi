extern "C"{
    void dgemv_cpp(int, double*, double*, double*);
    double ddot_cpp(int, double*, double*);
    void daxpy_cpp(int, double, double*, double*);
    double dnrm2_cpp(int, double*);
    void dscal_cpp(int, double, double*);
    void dhseqr_cpp(int, double*, double*, double*);
}

