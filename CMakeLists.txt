#
# Some help to write CMakeLists.cmake for fortran project:
#
# - http://www.cmake.org/Wiki/CMakeForFortranExample
# - http://www.ngssc.se/courses/specialized-courses/advanced-programming/make_cmake.pdf
# - https://github.com/SethMMorton/cmake_fortran_template
#

# require cmake 2.8.0
CMAKE_MINIMUM_REQUIRED(VERSION 2.8.0 FATAL_ERROR)

##################
project("laboetie")
##################

# Fortran is disabled by default, so we need to wake-up
# note: string given to enable_language is case-sensitive
enable_language(Fortran)

## guard against in-source builds
if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  message(FATAL_ERROR "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there. You may need to remove CMakeCache.txt. Follow instructions in README.md.")
endif()

# Build type :
# - RELEASE (with optimization flags, without openmp)
# - DEBUG   (with debug flags, without openmp)
# make sure that the default is RELEASE
if (NOT CMAKE_BUILD_TYPE)
    set (CMAKE_BUILD_TYPE RELEASE CACHE STRING
      "Choose the type of build, options are: None Debug Release."
      FORCE)
endif (NOT CMAKE_BUILD_TYPE)

# allow cmake to use custom modules located in cmake subdir
# (e.g. FindFFTW)
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
#find_package(MPI REQUIRED)

#
# check compiler (gfortran, ifort, ...)
# set FFLAGS depend on the compiler
# set where fortran module files will be located
#
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)
#SET (CMAKE_Fortran_COMPILER  gfortran)
message(STATUS "Using compiler ${Fortran_COMPILER_NAME}")

# Have the .mod files placed in the lib folder
SET(CMAKE_Fortran_MODULE_DIRECTORY "${PROJECT_BINARY_DIR}/mod")

if (Fortran_COMPILER_NAME MATCHES "gfortran.*")
  set (CMAKE_Fortran_FLAGS_RELEASE "-mcmodel=large -march=native -ffree-line-length-none -O3 -ffast-math -pedantic -fopenmp")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-mcmodel=large -march=native -ffree-line-length-none -O1             -pedantic -ffpe-trap=zero,underflow,overflow -Wall -fcheck=bounds -g")
  #-ftree-vectorizer-verbose=10 -fopt-info
  #set (CMAKE_Fortran_FLAGS_DEBUG   "-g -fbacktrace -pedantic -fwhole-file -Wline-truncation -Wcharacter-truncation -Wsurprising -Waliasing -fbounds-check -pg -frecursive -fcheck=all -Wall -ffpe-trap=zero,underflow,overflow")
elseif (Fortran_COMPILER_NAME MATCHES "ifort.*")
  set (CMAKE_Fortran_FLAGS_RELEASE "-O3")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
elseif (Fortran_COMPILER_NAME MATCHES "pgf90")
  set (CMAKE_Fortran_FLAGS_RELEASE "-O3")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
endif()

#############################
# Now build laboetie executable #
#############################

# set source files list
set(laboetie_SRCS
  src/module_propagation.f90
  src/module_bounceback.f90
  src/module_advect.f90
  src/backup_phi_c_plus_c_minus.f90
  src/charges_init.f90
  src/check_charge_conservation.f90
  src/module_external_forces.f90
  src/check_charge_distribution_equilibrium.f90
  src/module_tracers.f90
  src/electrostatic_pot.f90
  src/module_equilibration.f90
  src/init_simu.f90
  src/just_eq_smolu.f90
  src/main.f90
  src/module_collision.f90
  src/module_constants.f90
  src/module_geometry.f90
  src/module_input.f90
  src/module_io.f90
  src/module_lbmodel.f90
  src/module_mathematica.f90
  src/module_moment_propagation.f90
  src/module_myallocations.f90
  src/module_precision_kinds.f90
  src/module_system.f90
  src/module_time.f90
  src/poisson_nernst_planck.f90
  src/smolu.f90
  src/sor.f90
  src/supercell_definition.f90
#  src/module_mpi.F90
)



#
# Create a symlink to the input folder
#
EXECUTE_PROCESS(COMMAND "${CMAKE_COMMAND}" -E create_symlink "${CMAKE_SOURCE_DIR}/lb.in" "${CMAKE_BINARY_DIR}/lb.in")


if (UNIX)
  message(STATUS "Building on a UNIX-like platform")

  add_executable(laboetie ${laboetie_SRCS})

elseif(WIN32)
  message(STATUS "WIN32 platform not yet supported...")
endif()


#####################################
# Tell how to install this executable
#####################################

IF(WIN32)
    SET(CMAKE_INSTALL_PREFIX "C:\\Program Files")
ELSE()
    SET(CMAKE_INSTALL_PREFIX /usr/local)
ENDIF(WIN32)
INSTALL(TARGETS ${FOOEXE} RUNTIME DESTINATION bin)
