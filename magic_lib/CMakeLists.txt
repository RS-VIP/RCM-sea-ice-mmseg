cmake_minimum_required(VERSION 3.15)

project(magic_py VERSION 0.1 DESCRIPTION "MAGIC LIBRARY")

# ADD PYBIND11
# set(PYBIND11_PYTHON_VERSION 3.8 CACHE STRING "")
add_subdirectory(${CMAKE_SOURCE_DIR}/thirdparty/pybind11)
# set(PYBIND11_PYTHON_VERSION 3.8)

# modified permutohedral lattice gauss transform
#add_subdirectory(${CMAKE_SOURCE_DIR}/thirdparty/permuto)

# magic_core: contains c++ code for region-based image segmentation
add_subdirectory(${CMAKE_SOURCE_DIR}/magic_core)

# add the magic_py python module
pybind11_add_module(magic_py
        ${CMAKE_SOURCE_DIR}/magic_wrapper/src/PyRAG.cpp
#        ${CMAKE_SOURCE_DIR}/magic_wrapper/src/PyCurvatureGraph.cpp
        ${CMAKE_SOURCE_DIR}/magic_wrapper/src/PyGMM.cpp
        ${CMAKE_SOURCE_DIR}/magic_wrapper/src/Fisher.cpp
#        ${CMAKE_SOURCE_DIR}/magic_wrapper/src/permuto_wrapper.cpp
        ${CMAKE_SOURCE_DIR}/magic_wrapper/src/magic_py.cpp
        ${CMAKE_SOURCE_DIR}/magic_wrapper/src/TextureWrapper.cpp
        )

#target_link_libraries(magic_py PRIVATE magic_lib PERMUTO)
target_link_libraries(magic_py PRIVATE magic_lib)
target_include_directories(magic_py PUBLIC magic_wrapper/include)