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

cmake_minimum_required(VERSION 3.13)

project(rayleigh LANGUAGES C Fortran)

# load in version info and export it
file(STRINGS "${CMAKE_SOURCE_DIR}/VERSION" RAYLEIGH_PACKAGE_VERSION LIMIT_COUNT 1)


message(STATUS "")
message(STATUS "====================================================")
message(STATUS "===== Configuring Rayleigh ${RAYLEIGH_PACKAGE_VERSION}")
message(STATUS "====================================================")
message(STATUS "")

# Set up CMAKE_BUILD_TYPE to Debug if it is not set by the user.
if ("${CMAKE_BUILD_TYPE}" STREQUAL "")
  set(CMAKE_BUILD_TYPE
      "Release"
      CACHE STRING
      "Choose the type of build, options are: Debug, Release."
      FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
  "Debug" "Release")
endif()

message(STATUS "Configuring in \"" ${CMAKE_BUILD_TYPE} "\" mode.")
message(STATUS "")

# Set the default install directory, if run for the first time
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set_property(CACHE CMAKE_INSTALL_PREFIX PROPERTY VALUE "${CMAKE_SOURCE_DIR}")
endif()

message(STATUS "Installing executable in ${CMAKE_INSTALL_PREFIX}/bin")
message(STATUS "")

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

message(STATUS "===== Configuring external libraries ===============")
message(STATUS "")

find_package(MPI COMPONENTS C REQUIRED Fortran REQUIRED)
find_package(LAPACK COMPONENTS Fortran REQUIRED)
find_package(FFTW REQUIRED)

add_subdirectory(src)
