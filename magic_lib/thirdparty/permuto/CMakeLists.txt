cmake_minimum_required(VERSION 3.5)

project(PERMUTO VERSION 0.0.1 DESCRIPTION "modified permutohedral lattice stuff")

add_library(PERMUTO STATIC
        src/permutohedral.cpp
        )

set_target_properties(PERMUTO PROPERTIES
        VERSION ${PROJECT_VERSION}
        SOVERSION 0
        )

# no external dependencies
target_include_directories(PERMUTO PUBLIC include)
target_link_libraries(PERMUTO PUBLIC)