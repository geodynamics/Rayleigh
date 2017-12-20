# - Find FFTW3
# Find the native Fortran FFTW3 includes and library
#
#  FFTW_INCLUDES    - where to find fftw3.h
#  FFTW_LIBRARIES   - List of libraries when using FFTW.
#  FFTW_FOUND       - True if FFTW found.

if (FFTW3_Fortran_INCLUDES)
  # Already in cache, be silent
  set (FFTW3_Fortran_FIND_QUIETLY TRUE)
endif (FFTW3_Fortran_INCLUDES)

find_path (FFTW3_Fortran_INCLUDES fftw3.f)

find_library (FFTW3_LIBRARIES NAMES fftw3)

# handle the QUIETLY and REQUIRED arguments and set FFTW_FOUND to TRUE if
# all listed variables are TRUE
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (FFTW DEFAULT_MSG FFTW3_LIBRARIES FFTW3_Fortran_INCLUDES)

mark_as_advanced (FFTW3_LIBRARIES FFTW3_Fortran_INCLUDES)
