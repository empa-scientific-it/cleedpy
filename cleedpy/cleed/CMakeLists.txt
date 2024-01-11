cmake_minimum_required(VERSION 3.5)
project(cleed VERSION 1.0 LANGUAGES C)

# Specify the required C standard
set(CMAKE_C_STANDARD 90)
set(CMAKE_C_STANDARD_REQUIRED True)

# Set the output directories
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Where to find CMake modules
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

# Workaround for Apple Clang
# Disables the warning "implicit declaration of function '...' is invalid in C99"
if(APPLE AND "${CMAKE_C_COMPILER_ID}" STREQUAL "AppleClang")
    add_compile_options(-Wno-implicit-function-declaration)
endif()

# Set compiler and flags
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -O3 -march=native")

# Only set the compiler if strictly required; by default, let CMake find the most suitable one
#set(CMAKE_C_COMPILER "gcc-13")

# Global include directories
include_directories(include)

# Add the src directory
add_subdirectory(src)

# Build the test
add_executable(test_cleed test_cleed.c)

# Link the test
target_link_libraries(test_cleed cleed)
