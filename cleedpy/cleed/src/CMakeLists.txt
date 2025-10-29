# Add the libmat subdirectory
add_subdirectory(libmat)

# Compile cleed source files
file(GLOB CLEED_SRC "*.c")

# Remove a few src files that contain duplicated symbols
list(REMOVE_ITEM CLEED_SRC
    "${CMAKE_CURRENT_SOURCE_DIR}/linpphase.c"
    "${CMAKE_CURRENT_SOURCE_DIR}/lpcmktlnd_read.c"
)

add_library(cleed SHARED ${CLEED_SRC})

# Link libmat
target_link_libraries(cleed mat m dl)

# This is needed to tell the linker where to find `libmat`
# `libcleed` depends on it and setting RPATH allows us
# to resolve the correct path at runtime
# @loader_path and $ORIGIN point to the folder of the currently loaded shared library
if(APPLE)
    set_target_properties(cleed PROPERTIES
        INSTALL_RPATH "@loader_path"
        BUILD_WITH_INSTALL_RPATH TRUE
    )
elseif(UNIX)
    set_target_properties(cleed PROPERTIES
        INSTALL_RPATH "$ORIGIN"
        BUILD_WITH_INSTALL_RPATH TRUE
    )
endif()

# Install the cleed library into the Python package
# Use SKBUILD_PROJECT_NAME variable provided by scikit-build-core
install(TARGETS cleed
    LIBRARY DESTINATION ${SKBUILD_PROJECT_NAME}
    RUNTIME DESTINATION ${SKBUILD_PROJECT_NAME}
    COMPONENT cleed
)
