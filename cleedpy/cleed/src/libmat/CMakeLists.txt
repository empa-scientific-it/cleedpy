# Set BLAS/LAPACK vendor
set(BLA_VENDOR OpenBLAS)

# Automatically collect all source files matching the pattern "matXXX.c"
file(GLOB LIBMAT_SOURCES "mat*.c")

# Define the library from the collected source files
add_library(mat SHARED ${LIBMAT_SOURCES})

# Find BLAS
find_package(OpenBLAS REQUIRED)
include_directories(${OpenBLAS_INCLUDE_DIR})

# Link libmat with BLAS
target_link_libraries(mat ${OpenBLAS_LIB})
