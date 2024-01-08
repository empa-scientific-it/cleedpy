# Add the libmat subdirectory
add_subdirectory(libmat)

# Compile cleed source files
file(GLOB CLEED_SRC "*.c")
add_library(cleed ${CLEED_SRC})

# Link libmat
target_link_libraries(cleed libmat m dl)