# Set BLAS/LAPACK vendor (configurable via BLA_VENDOR CMake variable or environment)
# Supported vendors:
#   OpenBLAS,
#   Intel10_64lp (MKL),
#   FlexiBLAS,
#   ATLAS,
#   Generic,
#   All
# BLA_VENDOR is passed from pyproject.toml or can be set via: cmake -DBLA_VENDOR=OpenBLAS
if(NOT DEFINED BLA_VENDOR)
    set(BLA_VENDOR "All")
endif()

# Control static/dynamic linking (configurable via BLA_STATIC)
if(NOT DEFINED BLA_STATIC)
    set(BLA_STATIC OFF)
endif()

# Automatically collect all source files matching the pattern "matXXX.c"
file(GLOB LIBMAT_SOURCES "mat*.c")

# Define the library from the collected source files
add_library(mat SHARED ${LIBMAT_SOURCES})

# Find BLAS using CMake's built-in FindBLAS module
# This works with OpenBLAS, MKL, FlexiBLAS, ATLAS, Accelerate (macOS), and others
find_package(BLAS REQUIRED)

# FindBLAS does NOT provide include directories, so we need to find them manually
# Search for cblas.h in common locations
find_path(CBLAS_INCLUDE_DIR cblas.h
    PATHS
        /usr/include
        /usr/local/include
        /opt/homebrew/opt/openblas/include  # Homebrew ARM64 macOS
        /usr/local/opt/openblas/include     # Homebrew x86_64 macOS
        /opt/openblas/include
        /usr/include/openblas
        /usr/include/atlas
        $ENV{BLAS_HOME}/include
        $ENV{OpenBLAS_HOME}/include
        $ENV{MKLROOT}/include
    PATH_SUFFIXES
        openblas
        atlas
    DOC "Path to cblas.h"
)

# On macOS with Accelerate framework, cblas.h is in the SDK
if(APPLE AND NOT CBLAS_INCLUDE_DIR)
    # Accelerate framework headers are typically found by the compiler automatically
    set(CBLAS_INCLUDE_DIR "")
    message(STATUS "Using macOS Accelerate framework (headers found automatically)")
endif()

if(CBLAS_INCLUDE_DIR)
    message(STATUS "Found cblas.h in: ${CBLAS_INCLUDE_DIR}")
    target_include_directories(mat PRIVATE ${CBLAS_INCLUDE_DIR})
elseif(NOT APPLE)
    message(WARNING "Could not find cblas.h. Build may fail.")
endif()

# Link BLAS libraries
target_link_libraries(mat ${BLAS_LIBRARIES})

# Print BLAS information
message(STATUS "BLAS vendor: ${BLA_VENDOR}")
message(STATUS "BLAS libraries: ${BLAS_LIBRARIES}")

# Install the mat library into the Python package
install(TARGETS mat
    LIBRARY DESTINATION ${SKBUILD_PROJECT_NAME}
    RUNTIME DESTINATION ${SKBUILD_PROJECT_NAME}
    COMPONENT cleed
)
