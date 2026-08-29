cmake_minimum_required(VERSION 3.12)
project(mpas_jigsaw C)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

find_package(MPI REQUIRED COMPONENTS C)
find_package(PnetCDF REQUIRED)

add_executable(mkgrid mkgrid.c)
target_link_libraries(mkgrid PRIVATE PnetCDF::PnetCDF MPI::MPI_C m)

add_executable(write_hfun write_hfun.c)

install(TARGETS mkgrid write_hfun DESTINATION bin)
