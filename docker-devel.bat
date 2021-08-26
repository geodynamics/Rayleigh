:: echo off

docker run -it --rm -v $HOME:/work -e NOUIDWARN=1 geodynamics/rayleigh-devel-bionic:latest
