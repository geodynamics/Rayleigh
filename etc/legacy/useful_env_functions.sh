function gen_env_dirs {
    # "Returns" an array of all directories found in the user's environment
    #  The array is echo'd, and needs to be received into an array via syntax
    #  like that shown on the usage line below.
    #  Usage :   envdirs=($( gen_env_dirs ))
    local envfile="myenv"
    env > $envfile


    #I am keeping this here for notes right now.
    #IFS sets the internal field separator within read
    #The "or" ( || ) statement in the commented line appears to be
    #there in case read returns exit status 0, but a line of zero length is read.
    #while IFS='' read -r line || [[ -n "$line" ]]; do
    # Note also, from Jonanthon Anderson, that the double [[ above essentially invokes
    # "test" as with a single [, but it allows additional, bash-related functionality

    local dcounter=0
    while IFS='' read -r line 
    do
        
        # We are interested in directory names, not environment variable
        # names.  Separate each line into two pieces based on the "=" sign
        # Notes:
        #       ${line//\=/ } will replace ALL (the two //'s) instances of = with a space ' '
        #       Technically, we only want the first "=", but if this is a colon-separated
        #       list of directories, then there shouldn't be a second = anyway...
        #       The ( ) converts the result into an array (based on whitespace)
        local PIECES=(${line//\=/ })

        # Retain information to the right of the "=" sign
        local DIRTMP=${PIECES[1]}

        # Check to see if the first character is a "/", if so, this is a directory or list of directories
        if [[ ${DIRTMP:0:1} == "/" ]]
        then
            # separate out any colons
            for dir in ${DIRTMP//:/ }
            do
                ENVDIRS[$dcounter]=$dir
                ((dcounter++))
            done
        fi
    done < "$envfile"
    rm $envfile
    echo ${ENVDIRS[*]}
}

function parse_dirs {
    # Search a white-space separate list of directories ($1)
    # for those containing check_string in their name ($2).
    # "Return" and array via echo of those directories 
    # example usage:
    # mydirs=(/software /software/intel /software/intel/mkl)
    # intel_dirs=($( parse_dirs "${mydirs[*]}" "intel" ))
    # echo ${intel_dirs[*]}
    # (will print "/software/intel /softare/intel/mkl")
    local indirs=($1)
    local check_string=$2

    local dcounter=0
    for dir in ${indirs[*]}
    do
        if [[ $dir == *"$check_string"* ]]
        then        
            PDIRS[$dcounter]=$dir    
        fi
        ((dcounter++))
    done
    echo ${PDIRS[*]}

}
function verify_mklroot {
  # This function verifies that required MKL files are found in the
  # directory $1
  #We do this by checking for the existence of the 6 files below
  local FAILURE=FALSE
  local DUMMY=TRUE
  local CHECK_INC=$1/include/mkl_service.fi
  local CHECK_FFTW=$1/include/fftw/fftw3.f
  local CHECK_GF=$1/lib/intel64/libmkl_gf_lp64.so
  local CHECK_SEQ=$1/lib/intel64/libmkl_sequential.so
  local CHECK_CORE=$1/lib/intel64/libmkl_core.so

  if [[ ! -e "$CHECK_INC" ]]
  then
    FAILURE=TRUE
  fi

  if [[ ! -e "$CHECK_FFTW" ]]
    then
    FAILURE=TRUE
  fi

  if [[ ! -e "$CHECK_GF" ]]
    then
    FAILURE=TRUE
  fi

  if [[ ! -e "$CHECK_SEQ" ]]
    then
    FAILURE=TRUE
  fi

  if [[ ! -e "$CHECK_CORE" ]]
    then
    FAILURE=TRUE
  fi

  if [[ $FAILURE == "TRUE" ]]
  then
    echo $1
  else
    echo "0"
  fi
}
function find_mkl_root {
    local indirs=($1)
    #collect all directories that have "mkl" in their name
    mkldirs=($( parse_dirs "${indirs[*]}" "mkl" ))    
    #echo ${mkldirs[*]}
    local rcounter=0
    local dirtwo=''
    local dsplit=''
    #Loop over the directories with mkl in their name.
    #Check the full directory, as well as everything up through the "mkl"
    #This way if mkl/lib/intel64 in the environment, the

    #local FOUND_MKL="FALSE"
    local MKL_LOC="NONE"
    for dir in ${mkldirs[*]}
    do
        dsplit=(${dir//"mkl"/ })
        dirtwo=${dsplit[0]}"mkl"
        #echo ""
        #echo $dir
        #echo $dirtwo
        dcheck=$(verify_mklroot "$dir")
        dchecktwo=$(verify_mklroot "$dirtwo")
        #echo $dcheck $dchecktwo
        if [[ "$dcheck" == "0" ]]
        then
            MKL_LOC=$dir
        fi
        if [[ "$dchecktwo" == "0" ]]
        then
            MKL_LOC=$dirtwo
        fi

        ((rcounter++))
    done    
    echo $MKL_LOC
}
