#!/bin/bash
function find_mklroot {
  if [ -z "$1" ]
  then
    echo "RAYLEIGH MKL_ROOT is not set: "
    #This is where we will search for MKL in the environment
  else
    #VERIFY THAT CERTAIN NEEDED FILES ARE HERE
    MKL_LOC=$1
  fi
  FAILURE=FALSE
  DUMMY=TRUE
  CHECK_INC=$1/include/mkl_service.fi
  CHECK_FFTW=$1/include/fftw/fftw3.f
  CHECK_GF=$1/lib/intel64/libmkl_gf_lp64.so
  CHECK_SEQ=$1/lib/intel64/libmkl_sequential.so
  CHECK_CORE=$1/lib/intel64/libmkl_core.so

  if [[ -e "$CHECK_INC" ]]
  then
    DUMMY=TRUE
  else
    echo "ERROR --  Required file not found: "$CHECK_INC
    FAILURE=TRUE
  fi

  if [[ -e "$CHECK_FFTW" ]]
    then
    DUMMY=TRUE
  else
    FAILURE=TRUE
    echo "ERROR --  Required file not found: "$CHECK_FFTW
  fi

  if [[ -e "$CHECK_GF" ]]
    then
    DUMMY=TRUE
  else
    FAILURE=TRUE
    echo "ERROR --  Required file not found: "$CHECK_GF
  fi

  if [[ -e "$CHECK_SEQ" ]]
    then
    DUMMY=TRUE
  else
    FAILURE=TRUE
    echo "ERROR --  Required file not found: "$CHECK_SEQ
  fi

  if [[ -e "$CHECK_CORE" ]]
    then
    DUMMY=TRUE
  else
    FAILURE=TRUE
    echo "ERROR --  Required file not found: "$CHECK_CORE
  fi

  if [[ $FAILURE == "FALSE" ]]
  then
    echo $MKL_LOC
  else
    echo "FAILED to detect complete MKL installation"
  fi
}

FVERSION=GNU

if [[ $FVERSION == "INTEL" ]]
then
  RA_MKL=$(find_mklroot $RAYLEIGH_MKL_ROOT)
  INC_FLAGS="-I\$(RA_MKL)/include -I\$(RA_MKL)/include/fftw"
  echo $RA_MKL
  echo $INC_FLAGS
  LIB_FLAGS="-lstdc++ -mkl"
elif [[ $FVERSION == "GNU" ]]
then
  ##NOTE TO SELF:   NEED -m64 flag for gfortran in addition -- tucked into inc_flags below
  RA_MKL=$(find_mklroot $RAYLEIGH_MKL_ROOT)
  INC_FLAGS="-m64 -I\$(RA_MKL)/include -I\$(RA_MKL)/include/fftw"
  echo $RA_MKL
  echo $INC_FLAGS
  LIB_FLAGS="-L\$(RA_MKL)/lib/intel64 -Wl,--no-as-needed -lmkl_gf_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl"
  echo $LIB_FLAGS
fi


