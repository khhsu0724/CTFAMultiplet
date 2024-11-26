# Minimum CMake version
cmake_minimum_required(VERSION 3.15)

# Project configuration
project(CTFAMultiplet LANGUAGES CXX)

# Compiler settings
if(POLICY CMP0167)
  cmake_policy(SET CMP0167 NEW)
endif()

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -ffast-math -march=native")

message("-- CXX: ${CMAKE_CXX_COMPILER}")
message("-- CXXFLAGS: ${CMAKE_CXX_FLAGS}")

#-------------------------------------------------------------------
# Check BOOST
#-------------------------------------------------------------------
find_package(Boost REQUIRED)
if (Boost_FOUND)
    message(STATUS "Boost found!")
    message(STATUS "Include directory: ${Boost_INCLUDE_DIRS}")
    message(STATUS "Libraries: ${Boost_LIBRARIES}")
else()
    message(FATAL_ERROR "Boost not found. Please install Boost and try again.")
endif()


#-------------------------------------------------------------------
# Check OpenMP
#-------------------------------------------------------------------
if(APPLE)
    if(CMAKE_C_COMPILER_ID MATCHES "Clang")
        set(OpenMP_C "${CMAKE_C_COMPILER}")
        set(OpenMP_C_FLAGS "-pedantic -Xclang -fopenmp -Wno-unused-command-line-argument")
        set(OpenMP_C_LIB_NAMES "libomp" "libgomp" "libiomp5")
        set(OpenMP_libomp_LIBRARY ${OpenMP_C_LIB_NAMES})
        set(OpenMP_libgomp_LIBRARY ${OpenMP_C_LIB_NAMES})
        set(OpenMP_libiomp5_LIBRARY ${OpenMP_C_LIB_NAMES})
    endif()
    if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        set(OpenMP_CXX "${CMAKE_CXX_COMPILER}")
        set(OpenMP_CXX_FLAGS "-pedantic -Xclang -fopenmp -Wno-unused-command-line-argument")
        set(OpenMP_CXX_LIB_NAMES "libomp" "libgomp" "libiomp5")
        set(OpenMP_libomp_LIBRARY ${OpenMP_CXX_LIB_NAMES})
        set(OpenMP_libgomp_LIBRARY ${OpenMP_CXX_LIB_NAMES})
        set(OpenMP_libiomp5_LIBRARY ${OpenMP_CXX_LIB_NAMES})
    endif()
endif()

find_package(OpenMP REQUIRED)

if (OPENMP_FOUND)
    # include_directories("${OPENMP_INCLUDES}")
    # link_directories("${OPENMP_CXX_LIBRARIES}")
    message(STATUS "OpenMP found!")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    message(STATUS "OMP Include directory: ${OpenMP_CXX_INCLUDE_DIRS}")
    message(STATUS "OMP libraries: ${OpenMP_CXX_LIBRARIES}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif(OPENMP_FOUND)



#-------------------------------------------------------------------
# Check ARPACK Library
#-------------------------------------------------------------------
set(ARPACK_ROOT "" CACHE PATH "Path to the ARPACK installation")
if (ARPACK_ROOT)
    set(ARPACK_INCLUDE_DIR "${ARPACK_ROOT}/include")
    set(ARPACK_LIB_DIR "${ARPACK_ROOT}/lib64")
    # Validate paths
    if (NOT EXISTS "${ARPACK_LIB_DIR}")
        ### If lib64 not found, try again
        set(ARPACK_LIB_DIR "${ARPACK_ROOT}/lib")
    endif()
    if (EXISTS "${ARPACK_INCLUDE_DIR}" AND EXISTS "${ARPACK_LIB_DIR}")
        message(STATUS "Using user-provided ARPACK_ROOT: ${ARPACK_ROOT}")
        include_directories(${ARPACK_INCLUDE_DIR})
        link_directories(${ARPACK_LIB_DIR})
        find_library(ARPACK_LIBRARY
            NAMES arpack
            HINTS ${ARPACK_LIB_DIR}
        )
        message(STATUS "Found ARPACK include directory: ${ARPACK_INCLUDE_DIR}")
        message(STATUS "Found ARPACK library: ${ARPACK_LIBRARY}")
    else()
        message(FATAL_ERROR "Invalid ARPACK_ROOT: ${ARPACK_ROOT}")
    endif()
else()
    find_library(ARPACK_LIBRARY
        NAMES arpack
        HINTS /usr/local/lib /opt/homebrew/lib /usr/lib /usr/local/Cellar/ /usr/local/opt/
        HINTS /usr/lib/x86_64-linux-gnu/ $ENV{LD_LIBRARY_PATH}
    )

    find_path(ARPACK_INCLUDE_DIR
        NAMES arpack.h
        HINTS /usr/local/include /opt/homebrew/include /usr/include /usr/local/include/arpack 
        HINTS /usr/include/arpack $ENV{CPATH} 
        PATH_SUFFIXES arpack-ng # This is used for sherlock
    )

    # Verify that ARPACK was found
    if (NOT ARPACK_INCLUDE_DIR OR NOT ARPACK_LIBRARY)
        message(FATAL_ERROR 
            "ARPACK library and/or include files not found. Please ensure they are installed.\n"
            "If ARPACK is installed but not detected, please set the ARPACK_ROOT variable."
        )
    else()
        message(STATUS "Found ARPACK include directory: ${ARPACK_INCLUDE_DIR}")
        message(STATUS "Found ARPACK library: ${ARPACK_LIBRARY}")
    endif()
endif()



#-------------------------------------------------------------------
# Check LAPACK/MKL
#-------------------------------------------------------------------
message(STATUS "ENV{MKLROOT}: $ENV{MKLROOT}")
set(MKL_INTERFACE_FULL intel_lp64)
find_package(MKL PATHS $ENV{MKLROOT})
if (MKL_FOUND)
    set(MKL_INCLUDE_DIRS "$ENV{MKLROOT}/include")
    message(STATUS "Intel MKL found!")
    message(STATUS "Include directory: ${MKL_INCLUDE_DIRS}")
    message(STATUS "Libraries: ${MKL_LIBRARIES}")
else()
    message(WARNING "Intel MKL not found. Falling back to LAPACK.")
    # set(USE_MKL OFF) # NOT sure if I should comment this out?
endif()

if (NOT MKL_FOUND)
    find_package(LAPACK REQUIRED)
    if (LAPACK_FOUND)
        message(STATUS "LAPACK found!")
        message(STATUS "Libraries: ${LAPACK_LIBRARIES}")
    else()
        message(FATAL_ERROR "LAPACK not found. Ensure LAPACK is installed.")
    endif()
endif()

#-------------------------------------------------------------------
# Define directories and files
#-------------------------------------------------------------------
set(SRC_DIR ${CMAKE_SOURCE_DIR}/src)
set(BUILD_DIR ${CMAKE_BINARY_DIR})
set(INCLUDE_DIR ${SRC_DIR})

# Find source and header files
file(GLOB SRC_FILES ${SRC_DIR}/*.cpp)
file(GLOB HEADER_FILES ${SRC_DIR}/*.hpp)

# Include directories
include_directories(${INCLUDE_DIR})

# SET MKLROOT just in case 
set(CMAKE_INSTALL_RPATH "$ENV{MKLROOT}/lib;")
set(CMAKE_BUILD_RPATH "$ENV{MKLROOT}/lib;")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

#-------------------------------------------------------------------
# Create the executable
#-------------------------------------------------------------------
if(APPLE)
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lomp")
endif()
add_executable(CTFAMultiplet ${SRC_FILES})
set_target_properties(CTFAMultiplet PROPERTIES OUTPUT_NAME ./main-test)

#-------------------------------------------------------------------
# Link Libraries
#-------------------------------------------------------------------
target_link_libraries(CTFAMultiplet 
    PRIVATE
    ${OpenMP_LIBRARIES}  # OpenMP libraries
    ${ARPACK_LIBRARY}    # ARPACK library
    ${LAPACK_LIBRARIES}  # LAPACK libraries (fallback)
)


#-------------------------------------------------------------------
# Include Directories
#-------------------------------------------------------------------
target_include_directories(CTFAMultiplet 
    PRIVATE
    ${Boost_INCLUDE_DIRS}  # Boost headers
    ${ARPACK_INCLUDE_DIR}  # ARPACK headers
)

if (MKL_FOUND)
    target_link_libraries(CTFAMultiplet 
        PRIVATE $<IF:${MKL_FOUND},$<LINK_ONLY:MKL::MKL>,> 
    )
    target_include_directories(CTFAMultiplet 
        PRIVATE $<IF:${MKL_FOUND},$<TARGET_PROPERTY:MKL::MKL,INTERFACE_INCLUDE_DIRECTORIES>,> 
    )
    target_compile_options(CTFAMultiplet
        PRIVATE $<IF:${MKL_FOUND},$<TARGET_PROPERTY:MKL::MKL,INTERFACE_COMPILE_OPTIONS>,> 
    )
endif()


# Clean up build files (CMake takes care of this with built-in clean)


