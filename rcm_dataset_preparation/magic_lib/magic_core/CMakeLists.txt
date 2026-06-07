cmake_minimum_required(VERSION 3.15)
project(magic_lib VERSION 0.1 DESCRIPTION "MAGIC SHARED LIBRARY")

IF (WIN32)
  set(tinyxml2_SHARED_LIBS OFF) # flag to build tinyxml2 as a shared library
ELSE()
  set(tinyxml2_SHARED_LIBS ON) # flag to build tinyxml2 as a shared library
ENDIF()


# ADD TINYXML2
add_subdirectory(thirdparty/tinyxml2)
add_subdirectory(thirdparty/shapelib)

add_library(magic_lib STATIC

        # Common stuff
        Common/src/EstimateUtil.cpp

        # Images + basic image processing stuff
        Image/src/Image.cpp
        Image/src/MonoImage.cpp
        Image/src/Filter.cpp
        Image/src/GaussFilter.cpp
        Image/src/GradientCalculator.cpp
        Image/src/Watershed.cpp

        # Region adjacency graph and associated components
        Graph/src/AVLTree.cpp
        Graph/src/Graph.cpp
        Graph/src/Boundary.cpp
        Graph/src/Region.cpp
        Graph/src/GaussianRegion.cpp
        Graph/src/RegionAdjacencyGraph.cpp
        Graph/src/PartialEnumerationGraph.cpp
        Graph/src/StrengthPenaltyBoundary.cpp
        Graph/src/GaussianStrengthPenaltyBoundary.cpp
        Graph/src/GIEPGraph.cpp
        Graph/src/RegionMerger.cpp
        Graph/src/EdgeWeightKey.cpp

        # Matrix things
        Matrix/src/Matrix.cpp
        Matrix/src/Vector.cpp
        Matrix/src/Eigen.cpp

        # Various Classification stuff, mostly historic and not required for my use case. But lots of the
        # files depend on one another so I'd rather leave the structure in place instead of picking out
        # the parts that I need.
        Classifier/src/AnnealingUtil.cpp
        Classifier/src/RegionClassifier.cpp
        Classifier/src/RegionClustering.cpp
        Classifier/src/RegionBasedClustering.cpp
        Classifier/src/PixelBasedClustering.cpp
        Classifier/src/GaussianRegionClustering.cpp
        Classifier/src/MarkovNetClassifier.cpp
        Classifier/src/Sampler.cpp
        Classifier/src/GaussianClustering.cpp

        # Utilities for dealing with Polygons
        Polygon/src/PolygonDrawUtil.cpp
        Polygon/src/PolygonProjectionUtil.cpp

        # Texture feature extraction
        Texture/src/TextureFeatExtract.cpp
        Texture/src/SignalProcessing.cpp
        Texture/src/matrix_gabor.cpp

        )


IF (WIN32)
  # set stuff for windows
ELSE()
  # set stuff for linux
  target_compile_options(magic_lib PRIVATE "-fPIC")
ENDIF()


target_link_libraries(magic_lib PRIVATE shapelib tinyxml2)
target_include_directories(magic_lib

        PUBLIC
        Common/include
        Image/include
        Graph/include
        Matrix/include
        Classifier/include
        Radarsat/include
        Polygon/include
        Texture/include

        PUBLIC
        thirdparty/tinyxml2
        thirdparty/shapelib/include

        )

