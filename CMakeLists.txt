cmake_minimum_required(VERSION 3.12)
project(mpas_jigsaw C)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

find_package(MPI REQUIRED COMPONENTS C)
find_package(PnetCDF REQUIRED)

add_executable(mkgrid mkgrid.c)
target_compile_options(mkgrid PRIVATE -O3 -g)
target_link_libraries(mkgrid PRIVATE PnetCDF::PnetCDF MPI::MPI_C m)
