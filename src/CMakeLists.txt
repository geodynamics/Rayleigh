file(GLOB SOURCES
  Data_Structures/*.F90
  Diagnostics/*.F90
  IO/*.F90
  Math_Layer/*.F90
  Parallel_Framework/*.F90
  Physics/*.F90
  Test_Suite/*.F90
  IO/cmkdir.c)
list(REMOVE_ITEM SOURCES ${PROJECT_SOURCE_DIR}/src/Physics/Linear_Terms_Cart.F90)

add_custom_target(run_param_header COMMAND env FC=${CMAKE_Fortran_COMPILER} ${PROJECT_SOURCE_DIR}/src/gen_header.sh Run_Param_Header.F ${CMAKE_BUILD_TYPE} "compile flags unavailable" $<TARGET_PROPERTY:rayleigh,LINK_OPTIONS>)

add_executable(rayleigh ${SOURCES})

add_dependencies(rayleigh run_param_header)

target_include_directories(rayleigh PRIVATE ${PROJECT_SOURCE_DIR}/src/Include ${CMAKE_CURRENT_BINARY_DIR})

if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  target_compile_definitions(rayleigh PRIVATE GNU_COMPILER)
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
  target_compile_definitions(rayleigh PRIVATE INTEL_COMPILER)
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "Cray")
  target_compile_definitions(rayleigh PRIVATE CRAY_COMPILER)
endif()

if(MPI_Fortran_HAVE_F08_MODULE)
  target_compile_definitions(rayleigh PRIVATE USE_MPI_F08_BINDINGS)
endif()

target_link_libraries(rayleigh MPI::MPI_Fortran)
target_link_libraries(rayleigh LAPACK::LAPACK)
target_link_libraries(rayleigh FFTW::fftw3)

install(TARGETS rayleigh DESTINATION ${CMAKE_SOURCE_DIR}/bin)
