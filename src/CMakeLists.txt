#
#  Copyright (C) 2024 by the authors of the RAYLEIGH code.
#
#  This file is part of RAYLEIGH.
#
#  RAYLEIGH is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3, or (at your option)
#  any later version.
#
#  RAYLEIGH is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with RAYLEIGH; see the file LICENSE.  If not see
#  <http://www.gnu.org/licenses/>.

file(GLOB SOURCES
  Data_Structures/*.F90
  Diagnostics/*.F90
  IO/*.F90
  Math_Layer/*.F90
  Parallel_Framework/*.F90
  Physics/*.F90
  Test_Suite/*.F90
  IO/cmkdir.c)
list(REMOVE_ITEM SOURCES ${PROJECT_SOURCE_DIR}/src/Physics/Linear_Terms_Cart.F90)

add_custom_target(run_param_header COMMAND env FC=${CMAKE_Fortran_COMPILER} ${PROJECT_SOURCE_DIR}/src/gen_header.sh Run_Param_Header.F $<$<STREQUAL:${CMAKE_BUILD_TYPE},>:Default>${CMAKE_BUILD_TYPE} FROM_CMAKE $<TARGET_PROPERTY:rayleigh,LINK_OPTIONS>)

add_executable(rayleigh ${SOURCES})

add_dependencies(rayleigh run_param_header)

target_include_directories(rayleigh PRIVATE ${PROJECT_SOURCE_DIR}/src/Include ${CMAKE_CURRENT_BINARY_DIR})

set(RAYLEIGH_CPU_OPTIMIZATIONS "none" CACHE STRING "Specify a CPU architecture to build for. Currently only 'native' or 'none' are supported.")
set_property(CACHE RAYLEIGH_CPU_OPTIMIZATIONS PROPERTY STRINGS "none" "native")

if (RAYLEIGH_CPU_OPTIMIZATIONS STREQUAL "native")
  if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    target_compile_options(rayleigh PRIVATE -march=native)
    message(STATUS "enabling native optimizations")
  elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
    target_compile_options(rayleigh PRIVATE -xHost)
  elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "IntelLLVM")
    target_compile_options(rayleigh PRIVATE -xHost)
  else()
    message(FATAL_ERROR "For your compiler no native optimizations are automatically detected. Please specify compiler flags manually and set 'CPU_OPTIMIZATIONS' to 'none'")
  endif()
elseif (RAYLEIGH_CPU_OPTIMIZATIONS STREQUAL "none")
  # do nothing, respect manually set compiler flags
else()
  message(FATAL_ERROR "Only 'native' and 'none' are supported for parameter 'CPU_OPTIMIZATIONS'")
endif()

if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  target_compile_definitions(rayleigh PRIVATE GNU_COMPILER)
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
  target_compile_definitions(rayleigh PRIVATE INTEL_COMPILER)
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "Cray")
  target_compile_definitions(rayleigh PRIVATE CRAY_COMPILER)
endif()

if(MPI_Fortran_HAVE_F08_MODULE)
  target_compile_definitions(rayleigh PRIVATE USE_MPI_F08_BINDINGS)
endif()

target_link_libraries(rayleigh MPI::MPI_Fortran)
target_link_libraries(rayleigh LAPACK::LAPACK)
target_link_libraries(rayleigh FFTW::fftw3)

message(STATUS "")

install(TARGETS rayleigh RUNTIME DESTINATION bin COMPONENT runtime)
