cmake_minimum_required(VERSION 3.10)
project(MatrixProject)

set(CMAKE_CXX_STANDARD 17)

add_executable(matrix_program main.cpp)  # Убираем Matrix.cpp

if (MSVC)
    target_compile_options(matrix_program PRIVATE /W4)
else()
    target_compile_options(matrix_program PRIVATE -Wall -Wextra -Wpedantic)
endif()
