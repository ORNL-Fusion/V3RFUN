cmake_minimum_required (VERSION 3.14)

project (v3rfun C Fortran)

find_package (OpenMP)

add_executable (xv3rfun)
add_dependencies (xv3rfun stell)

target_link_libraries (xv3rfun stell)
target_include_directories (xv3rfun PUBLIC $<TARGET_PROPERTY:stell,BINARY_DIR>)

add_subdirectory (Sources)

if (OpenMP_Fortran_FOUND)
    target_link_libraries (xv3rfun OpenMP::OpenMP_Fortran)
endif ()
