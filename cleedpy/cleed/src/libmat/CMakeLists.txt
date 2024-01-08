# Set BLAS/LAPACK vendor
set(BLA_VENDOR OpenBLAS)

# Automatically collect all source files matching the pattern "matXXX.c"
file(GLOB LIBMAT_SOURCES "mat*.c")

# Define the library from the collected source files
add_library(libmat ${LIBMAT_SOURCES})

# Include directories for BLAS/LAPACK
include_directories(/usr/include)

# Find BLAS/LAPACK
find_package(BLAS REQUIRED)
#find_package(LAPACK REQUIRED)

# Link libmat with BLAS and LAPACK
target_link_libraries(libmat ${BLAS_LIBRARIES})

