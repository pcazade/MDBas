cmake_minimum_required (VERSION 2.6)

set(TGT "MDBas" CACHE STRING "The executable name.")

project(${TGT} C)

# check if MPI is required by the user
option(ENABLE_MPI  "Enable MPI Parallel version" OFF)
if(ENABLE_MPI)
    find_package(MPI REQUIRED)
endif(ENABLE_MPI)

# list of include folders with .h files
include_directories(include dSFMT)

# list all source files which are common to all releases
set(
COMMON_SRCS
src/elec.c
src/errors.c
src/init.c
src/internal.c
src/list.c
src/memory.c
src/numderiv.c
src/spme.c
src/user.c
src/vdw.c
src/energy.c
src/ewald.c
src/integrate.c
src/io.c 
src/mdbas.c
src/minim.c
src/rand.c
src/shake.c
src/timing.c
src/utils.c
dSFMT/dSFMT.c
)

# list sources only used for a serial CPU execution
set(SERIAL_SRCS  src/serial.c)

# list sources only used fir a parallel execution using MPI
set(PARALLEL_SRCS src/parallel.c)

# BROKEN : OpenCL part
# set( 
# src/ocl/initOCL.c 
# src/ocl/prepareOCL.c 
# src/ocl/vdwOCL.c.BROKEN
# ${COMMON_SRCS} OPENCL_SRCS)

# add custom options for compiler, but common to all (serial or parallel ... ) targets.
set(CUSTOM_C_FLAGS "-std=gnu99 -Wall -Wextra " CACHE STRING "Custom C compiler options.")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${CUSTOM_C_FLAGS}")

#common defines
add_definitions(-DHAVE_SSE2 -DDSFMT_MEXP=19937 -DTIMING -DFFTW)

# common libraries
link_libraries(m rt dl fftw3)

# check if we need to add MPI related files
if(NOT ENABLE_MPI)

    set(SRCS ${COMMON_SRCS} ${SERIAL_SRCS})
    add_executable(${TGT} ${SRCS})
    
else(NOT ENABLE_MPI)

    set(SRCS ${COMMON_SRCS} ${PARALLEL_SRCS})
    include_directories(${MPI_INCLUDE_PATH})
    add_executable(${TGT}_MPI ${SRCS})
    target_link_libraries(${TGT}_MPI ${MPI_LIBRARY})
    add_definitions(-DUSING_MPI)
    
endif(NOT ENABLE_MPI)


