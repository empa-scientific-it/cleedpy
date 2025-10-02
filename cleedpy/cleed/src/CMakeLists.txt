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

# Install the cleed library into the Python package
# scikit-build-core will place this in the wheel at cleedpy/cleed/lib/
install(TARGETS cleed
    LIBRARY DESTINATION cleed/lib
    RUNTIME DESTINATION cleed/bin
    COMPONENT cleed
)
