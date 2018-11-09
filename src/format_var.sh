#!/bin/bash
#########################################################################
#  format_var.sh
#
#  For use in conjunction with gen_header.sh.
#  Expands a lengthy string into multiple lines of 
#  length 55 or less.  Appends each line with appropriate
#  Fortran line-continuation and string-concatenation
#  characters.
#
#  Calling syntax:
#  ./format_var.sh $VNAME $OFILE
#
#
#  Command-line inputs:
#       VNAME :  The variable declaration to expand into multi-line format 
#       OFILE :  The .F that each line is written to
###########################################################################

UNAME=$1
OFILE=$2

vlen=${#UNAME}

COUNTER=0
SLEN=55

i=0
IFS='%'
while [ $i -lt $vlen ]; do
    nextline=${UNAME:$i:$SLEN}

    if [ $i -eq 0 ]
    then
        nextline="$nextline"\"
    else
        nextline="        "\""$nextline"\"
    fi

    if [ $(( $i + $SLEN )) -ge $vlen ]
    then
        nextline="$nextline"\/\/"\"---\""\/\/"Char(0)"
    else
        nextline="$nextline"\/\/" &"
    fi
    echo $nextline >> $OFILE
    let i=i+$SLEN
done
echo '' >> $OFILE
unset IFS
