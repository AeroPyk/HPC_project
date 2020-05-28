cmake_minimum_required(VERSION 3.16)
project(ProjectMM C)

set(CMAKE_C_FLAGS -fopenmp)
set(CMAKE_C_STANDARD 99)

find_package(MPI REQUIRED) # define MPI function

add_executable(main main.c ext.h ext.c)

target_link_libraries(main ${MPI_LIBRARIES})