cmake_minimum_required(VERSION 3.7)

project(shapelib VERSION 1.5.0 DESCRIPTION "Shapelib Library")

include_directories(include)

add_library(shapelib STATIC
        src/dbfopen.c
        src/safileio.c
        src/shpopen.c
        )

set_target_properties(shapelib PROPERTIES
        VERSION ${PROJECT_VERSION}
        SOVERSION 0
        )

# no external dependencies
target_include_directories(shapelib PUBLIC include)
target_link_libraries(shapelib PUBLIC)