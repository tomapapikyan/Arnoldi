add_executable(main main.cpp)

target_link_libraries(main CppBLAS VectorFunctions MatrixTemplate)

#find_package(OpenMP)

#target_link_libraries(main PRIVATE OpenMP::OpenMP_CXX)
