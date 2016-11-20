# Minimum cmake version
cmake_minimum_required (VERSION 3.5.0)

# Name of project and that it is C++ only.
set(code_name "tfo_finder")
project (${code_name} )

# ----------------------------------------------------------------------------
# Dependencies
# ---------------------------------------------------------------------------- 
# Search for zlib as a dependency for SeqAn.
find_package (ZLIB)
find_package (Boost REQUIRED COMPONENTS regex)
if (NOT Boost_FOUND )
        message (STATUS "Boost not found, not building.")
            return ()
        endif ()
# Load the SeqAn module and fail if not found.
# Would be nice to add an ifdef depending on the machie
execute_process(COMMAND hostname OUTPUT_VARIABLE computer_name)
if (computer_name MATCHES "^air")
    message (STATUS "We are in air, configuring accordingly" )
    set (Vienna_LIBRARIES "/usr/loca/lib/")
    include_directories ("/usr/local/include/ViennaRNA" )
    message ( STATUS "This is the compiler ${CMAKE_CXX_COMPILER}" )
    message ( STATUS "Make sure to use -D CMAKE_CXX_COMPILER is set to 3.9 in air " )
    message ( STATUS "This is the compiler ${CMAKE_C_COMPILER}" )
    message ( STATUS "Make sure to use -D CMAKE_C_COMPILER is set to 3.9 in air " )
else() 
# assume you are in the iMac
    message (STATUS "We might be in iMac, configuring accordingly" )
    set (Vienna_LIBRARIES "/Users/guillem/work/soft_local_copy/ViennaRNA/include/ViennaRNA/")
    include_directories ("/Users/guillem/work/soft_local_copy/ViennaRNA/include/ViennaRNA" )
    message ( STATUS "This is the compiler ${CMAKE_CXX_COMPILER}" )
    message ( STATUS "Make sure to use -D CMAKE_CXX_COMPILER is set to 3.8 in air " )
    message ( STATUS "This is the compiler ${CMAKE_C_COMPILER}" )
    message ( STATUS "Make sure to use -D CMAKE_C_COMPILER is set to 3.8 in air " )
    set(CMAKE_C_COMPILER /opt/local/bin/clang-mp-3.8 )
    set(CMAKE_CXX_COMPILER /opt/local/bin/clang++-mp-3.8 )
endif()
find_package (SeqAn REQUIRED)
set (OpenMP_C_FLAGS /opt/local/lib/libomp )
set (OpenMP_CXX_FLAGS /opt/local/lib/libomp )
find_package (OpenMP REQUIRED)


# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

#set ( EIGEN3_INCLUDE_DIR "/opt/local/include/eigen3/Eigen" )
#if ( NOT EIGEN3_INCLUDE_DIR )
#    message( FATAL_ERROR "Please point the environment variable EIGEN3_INCLUDE_DIR to the include directory of your Eigen3 installation.")
#endif()
# Add include directories.
include_directories (${SEQAN_INCLUDE_DIRS} ${EIGEN3_INCLUDE_DIR} )

# Add definitions set by find_package (SeqAn).
add_definitions (${SEQAN_DEFINITIONS})

# Add CXX flags found by find_package (SeqAn).

# Do use -O3, makes the code way faster
set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -I/opt/local/include/libomp -fopenmp ")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS} -g -I/opt/local/include/libomp  -fopenmp -O3  ")

# Add executable and link against SeqAn dependencies.
add_executable (${code_name} ${code_name}.cpp nn_tm.hpp constants_nn.hpp commandline_parse.hpp find_tfo.hpp do_rnafold.c )
#add_library ( do_rnafold do_rnafold.c do_rnafold.h) 
find_library (RESULT libRNA.a PATHS ${Vienna_LIBRARIES} )
message (STATUS ${RESULT} )
target_link_libraries (${code_name} ${SEQAN_LIBRARIES} ${OpenPM_LIBRARIES} ${Boost_LIBRARIES} RNA )

