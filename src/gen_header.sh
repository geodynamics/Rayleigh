#!/bin/bash
#########################################################################
#  gen_header.sh
#
#  Helper utility used for generating the Run_Parameter include file
#  during Rayleigh's configuration process.
#
#  Creates a header file containing string variable declarations initialized
#  to values describing the build details.
#
#  Calling syntax:
#  ./gen_header.sh $HFILE "$RVER" "$FFLAGS" "$LIBS"
#
#
#  Command-line inputs:
#       HFILE  :  The include .F to create
#	RVER   :  The version of Rayleigh being compiled (e.g., opt, dbg)
#	FFLAGS :  The compile flags for this version
#	LIBS   :  The library-linking flags used to comipile this version.
###########################################################################
HFILE=$1


rm -rf $HFILE

#Git Hash
pref="    Character(len=128) :: build_git_commit = \"---git hash---"
val=$(git log -1 --pretty=format:"%H +++ %cd")
var="$pref$val"

./src/format_var.sh "$var" $HFILE

#Git url
pref="    Character(len=128) :: build_git_url = \"---git url---"
val=$(git config remote.origin.url)
var="$pref$val"
./src/format_var.sh "$var" $HFILE

#Git branch
pref="    Character(len=128) :: build_git_branch = \"---git branch---"
val=$(git rev-parse --abbrev-ref HEAD)
var="$pref$val"
./src/format_var.sh "$var" $HFILE

#Build date
pref="    Character(len=128) :: build_date = \"---build date---"
val=$(date)
var="$pref$val"
./src/format_var.sh "$var" $HFILE

#Machine name
pref="    Character(len=256) :: build_machine = \"---build machine---"
val=$(uname -a)
var="$pref$val"
./src/format_var.sh "$var" $HFILE

#Build directory
pref="    Character(len=512) :: build_dir = \"---build dir---"
val=$PWD
var="$pref$val"
./src/format_var.sh "$var" $HFILE

#Custom directory specified (if any)
pref="    Character(len=512) :: build_custom_dir = \"---custom dir---"
val=$CUSTOMROOT
var="$pref$val"
./src/format_var.sh "$var" $HFILE

#Build version
pref="    Character(len=64 ) :: build_version = \"---build version---"
val=$2
var="$pref$val"
./src/format_var.sh "$var" $HFILE

#Build flags
pref="    Character(len=512) :: build_fflags = \"---build fflags---"
val="$3"
var="$pref$val"
./src/format_var.sh "$var" $HFILE

#Build libraries
pref="    Character(len=512) :: build_lib = \"---build lib---"
val="$4"
var="$pref$val"
./src/format_var.sh "$var" $HFILE

#Fortran compiler version
pref="    Character(len=512) :: FC_version = \"---FC version---"
ver=$($FC --version)
val=$(echo "$ver" | sed -n 1p) # exclude lines 2 onward (copyright info)
var="$pref$val"
./src/format_var.sh "$var" $HFILE

#Fortran compiler location
pref="    Character(len=512) :: FC_location = \"---FC location---"
val=$(which $FC)
var="$pref$val"
./src/format_var.sh "$var" $HFILE

#more $HFILE
