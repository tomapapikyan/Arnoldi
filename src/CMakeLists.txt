add_library(Arnoldi ../include/Arnoldi_Iterations.h Arnoldi_Iterations.cpp)

add_library(CppBLAS ../include/CppBLAS.h CppBLAS.f90)

add_library(MatVec ../include/MatVec.h MatVec.cpp)

find_package(BLAS REQUIRED)

target_link_libraries(CppBLAS ${BLAS_LIBRARIES})

find_package(LAPACK REQUIRED)

target_link_libraries(CppBLAS ${LAPACK_LIBRARIES})

find_package(OpenMP)

target_link_libraries(Arnoldi PRIVATE OpenMP::OpenMP_CXX)
