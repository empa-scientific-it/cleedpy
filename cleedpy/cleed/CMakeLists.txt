cmake_minimum_required(VERSION 3.0)
project(cleed VERSION 1.0 LANGUAGES C)

# Specify the required C standard
set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED True)

# Global include directories or compiler options here
include_directories(include)

# Add the src directory
add_subdirectory(src)

# Build the test
add_executable(test_cleed test_cleed.c)

# Link
target_link_libraries(test_cleed cleed)