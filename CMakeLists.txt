cmake_minimum_required (VERSION 3.14)

project (v3rfun C Fortran)

find_package (OpenMP)

add_executable (xv3rfun)

target_link_libraries (xv3rfun 

                       PUBLIC 

                       stell
                       $<$<BOOL:xv3rfun>:OpenMP::OpenMP_Fortran>
)

add_subdirectory (Sources)
