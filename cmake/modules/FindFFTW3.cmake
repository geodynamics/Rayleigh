# copied from here: https://github.com/jedbrown/cmake-modules/blob/master/FindFFTW.cmake


# - Find FFTW
# Find the native FFTW includes and library
#
#  FFTW_INCLUDES    - where to find fftw3.h
#  FFTW_LIBRARIES   - List of libraries when using FFTW.
#  FFTW_FOUND       - True if FFTW found.

IF (FFTW_INCLUDES)
  # Already in cache, be silent
  SET (FFTW_FIND_QUIETLY TRUE)
ENDIF (FFTW_INCLUDES)

FIND_PATH (FFTW_INCLUDES fftw3.h)
INCLUDE_DIRECTORIES(${FFTW_INCLUDES})

FIND_LIBRARY (FFTW_LIBRARIES NAMES fftw3)

# handle the QUIETLY and REQUIRED arguments and set FFTW_FOUND to TRUE if
# all listed variables are TRUE
INCLUDE (FindPackageHandleStandardArgs)
find_package_handle_standard_args (FFTW DEFAULT_MSG FFTW_LIBRARIES FFTW_INCLUDES)

MARK_AS_ADVANCED (FFTW_LIBRARIES FFTW_INCLUDES)
