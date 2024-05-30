#add_library(Arnoldi ../include/Arnoldi_Iterations.h Arnoldi_Iterations.cpp)

add_library(CppBLAS ../include/CppBLAS.h CppBLAS.f90)

add_library(VectorFunctions ../include/VectorFunctions.h VectorFunctions.cpp)

add_library(MatrixTemplate ../include/MatrixTemplate.h MatrixTemplate.cpp)

find_package(LAPACK REQUIRED)

target_link_libraries(CppBLAS ${LAPACK_LIBRARIES})

find_package(OpenMP)

#target_link_libraries(Arnoldi PRIVATE OpenMP::OpenMP_CXX)

target_link_libraries(VectorFunctions PRIVATE OpenMP::OpenMP_CXX)

target_link_libraries(MatrixTemplate PRIVATE OpenMP::OpenMP_CXX)
