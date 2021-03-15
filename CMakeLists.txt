cmake_minimum_required (VERSION 3.14)

project (makegrid C Fortran)

find_package (OpenMP)

add_executable (v3rfun)
add_dependencies (v3rfun stell)

target_link_libraries (v3rfun stell)
target_include_directories (v3rfun PUBLIC $<TARGET_PROPERTY:stell,BINARY_DIR>)

add_subdirectory (Sources)

if (OpenMP_Fortran_FOUND)
    target_link_libraries (v3rfun OpenMP::OpenMP_Fortran)
endif ()
