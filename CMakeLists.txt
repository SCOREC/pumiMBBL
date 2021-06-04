cmake_minimum_required(VERSION 3.19)
project(pumiMBBL CXX )

find_package(Kokkos)

add_subdirectory (src)
add_subdirectory (test)
