# Set BLAS/LAPACK vendor
set(BLA_VENDOR OpenBLAS)

# Automatically collect all source files matching the pattern "matXXX.c"
file(GLOB LIBMAT_SOURCES "mat*.c")

# Define the library from the collected source files
add_library(libmat STATIC ${LIBMAT_SOURCES})

# Find BLAS/LAPACK
find_package(BLAS REQUIRED)

# Include directories for BLAS/LAPACK
find_file(BLAS_HEADER "cblas.h" HINTS ${CMAKE_PREFIX_PATH})
if(BLAS_HEADER)
    get_filename_component(BLAS_INCLUDE_DIR ${BLAS_HEADER} DIRECTORY)
    include_directories(${BLAS_INCLUDE_DIR})
endif()

# Link libmat with BLAS/LAPACK
target_link_libraries(libmat ${BLAS_LIBRARIES})
