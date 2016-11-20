# Minimum cmake version
cmake_minimum_required (VERSION 3.0.0)

# Name of project and that it is C++ only.
set(code_name "approxseqmatch")
project (${code_name} CXX)

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

# Search for zlib as a dependency for SeqAn.
find_package (ZLIB)
find_package (Boost REQUIRED COMPONENTS iostreams-mt)
# Stop here if we cannot find Boost or OpenMP.
if (NOT Boost_FOUND )
        message (STATUS "  Boost not found, not building fiona.")
            return ()
        endif ()
# Would be nice to add an ifdef depending on the machie
execute_process(COMMAND hostname OUTPUT_VARIABLE computer_name)
if (computer_name MATCHES "^air")
    message (STATUS "We are in air, configuring accordingly" )
    message ( STATUS "This is the compiler ${CMAKE_CXX_COMPILER}" )
    message ( STATUS "Make sure to use -D CMAKE_CXX_COMPILER is set to 3.9 in air " )
    message ( STATUS "This is the compiler ${CMAKE_C_COMPILER}" )
    message ( STATUS "Make sure to use -D CMAKE_C_COMPILER is set to 3.9 in air " )
else()
# assume you are in the iMac
    message (STATUS "We might be in iMac, configuring accordingly" )
    message ( STATUS "This is the compiler ${CMAKE_CXX_COMPILER}" )
    message ( STATUS "Make sure to use -D CMAKE_CXX_COMPILER is set to 3.8 in air " )
    message ( STATUS "This is the compiler ${CMAKE_C_COMPILER}" )
    message ( STATUS "Make sure to use -D CMAKE_C_COMPILER is set to 3.8 in air " )
    set(CMAKE_C_COMPILER /opt/local/bin/clang-mp-3.8 )
    set(CMAKE_CXX_COMPILER /opt/local/bin/clang++-mp-3.8 )
endif()

# Load the SeqAn module and fail if not found.
find_package (SeqAn REQUIRED)

# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

# Add include directories.
include_directories (${SEQAN_INCLUDE_DIRS})
include_directories (SYSTEM ${Boost_INCLUDE_DIR})

# Add definitions set by find_package (SeqAn).
add_definitions (${SEQAN_DEFINITIONS})

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS} ${Boost_CXX_FLAGS}")

# Add executable and link against SeqAn dependencies.
add_executable (${code_name} ${code_name}.cpp)
target_link_libraries (${code_name} ${SEQAN_LIBRARIES} ${Boost_LIBRARIES})
