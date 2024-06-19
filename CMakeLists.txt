cmake_minimum_required(VERSION 3.13)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

project(rayleigh LANGUAGES C Fortran)

find_package(MPI COMPONENTS C REQUIRED Fortran REQUIRED)
find_package(LAPACK COMPONENTS Fortran REQUIRED)
find_package(FFTW REQUIRED)

add_subdirectory(src)
