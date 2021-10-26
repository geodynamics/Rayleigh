:: echo off

docker run -it --rm -v "%USERPROFILE%":/work -e NOUIDWARN=1 geodynamics/rayleigh-devel-bionic:latest
