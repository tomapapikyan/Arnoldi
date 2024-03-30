cmake_minimum_required(VERSION 3.20)

project (MAIN
         LANGUAGES CXX Fortran)

include_directories(include)

add_subdirectory(src)

add_subdirectory(test)
