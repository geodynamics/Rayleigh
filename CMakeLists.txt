# CMake project file for RAYLEIGH

##################################################
# Define the project and the depencies that it has
##################################################

CMAKE_MINIMUM_REQUIRED(VERSION 3.0.0)
PROJECT(RAYLEIGH Fortran C)

# Set the Rayleigh version
SET(VERSION 0.9)

# Add our local modules to the module path
SET(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/Modules/")

IF(NOT CMAKE_Fortran_COMPILER_SUPPORTS_F90)
  MESSAGE(FATAL_ERROR "Fortran compiler does not support F90")
ENDIF(NOT CMAKE_Fortran_COMPILER_SUPPORTS_F90)

IF (NOT MPI_Fortran_FOUND)
  FIND_PACKAGE (MPI REQUIRED)
ENDIF (NOT MPI_Fortran_FOUND)

IF ( MPI_FOUND )
   SET (CMAKE_Fortran_COMPILER  mpif90)
   SET (CMAKE_CC_COMPILER  mpicc)
   SET (CMAKE_CXX_COMPILER  mpicxx)
ENDIF()

set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake/modules)
FIND_PACKAGE(FFTW3)


############################################################
# Define the actual files and folders that make up the build
############################################################

# Define the executable name
SET(RayleighEXE "rayleigh.dbg")

# Define some directories
SET(SRCDIR ${CMAKE_SOURCE_DIR}/src)
SET(LIBDIR ${CMAKE_SOURCE_DIR}/lib)
SET(BINDIR ${CMAKE_SOURCE_DIR}/bin)

# Have the .mod files placed in the lib folder
SET(CMAKE_Fortran_MODULE_DIRECTORY ${LIBDIR})

# The source for the Rayleigh binary and have it placed in the bin folder
ADD_SUBDIRECTORY(src)

# Add a distclean target to the Makefile
ADD_CUSTOM_TARGET(distclean 
    COMMAND ${CMAKE_COMMAND} -P ${CMAKE_SOURCE_DIR}/distclean.cmake
)
