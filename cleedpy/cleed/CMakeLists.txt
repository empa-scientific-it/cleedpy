cmake_minimum_required(VERSION 3.5)
project(cleed VERSION 1.0 LANGUAGES C)

# Specify the required C standard
set(CMAKE_C_STANDARD 90)
set(CMAKE_C_STANDARD_REQUIRED True)

# Set the output directories
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

# Where to find CMake modules
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

# Workaround for Apple Clang
# Disables the warning "implicit declaration of function '...' is invalid in C99"
if(APPLE AND "${CMAKE_C_COMPILER_ID}" STREQUAL "AppleClang")
    add_compile_options(-Wno-implicit-function-declaration)
endif()

# Set compiler optimization flags
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -O3")

# Architecture-specific optimizations (optional)
# By default, build portable binaries. Set CLEEDPY_ARCH_NATIVE=ON for -march=native
# Or set CLEEDPY_ARCH_FLAGS to specify custom architecture flags
if(CLEEDPY_ARCH_NATIVE)
    message(STATUS "Building with -march=native (optimized for build machine)")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -march=native")
elseif(DEFINED CLEEDPY_ARCH_FLAGS AND NOT CLEEDPY_ARCH_FLAGS STREQUAL "")
    message(STATUS "Building with custom architecture flags: ${CLEEDPY_ARCH_FLAGS}")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${CLEEDPY_ARCH_FLAGS}")
else()
    # Portable build - no architecture-specific flags
    message(STATUS "Building portable binary (no architecture-specific optimizations)")
    message(STATUS "  - For AMD Rome (x86_64): set CLEEDPY_ARCH_FLAGS='-march=znver2' or '-march=x86-64-v3'")
    message(STATUS "  - For Grace-Hopper (aarch64): set CLEEDPY_ARCH_FLAGS='-mcpu=neoverse-v2' or '-march=armv8.2-a'")
    message(STATUS "  - For native build: set CLEEDPY_ARCH_NATIVE=ON")
endif()

# Only set the compiler if strictly required; by default, let CMake find the most suitable one
#set(CMAKE_C_COMPILER "gcc-13")

# Global include directories
include_directories(include)

# Add the src directory
add_subdirectory(src)
