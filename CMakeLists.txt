# Specify the minimum version for CMake

cmake_minimum_required(VERSION 3.7)
set(CMAKE_BUILD_TYPE Debug)

# Project's name

project(plotTools)

find_package(ROOT 6 REQUIRED COMPONENTS Minuit)
include(${ROOT_USE_FILE})

# Set the output folder where your program will be created
set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/bin)
set(LIBRARY_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/lib)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})

#EIGEN3 package needed if dont have one, go git checkout and put the package path here
set(EIGEN3DIR /Users/Liby/work/tools/eigen)

# The following folder will be included
FILE(GLOB INCLUDES ${PROJECT_SOURCE_DIR}/include/ ${PROJECT_SOURCE_DIR}/include/*)
include_directories(${INCLUDES} ${PROJECT_SOURCE_DIR}/include ${EIGEN3DIR})

AUX_SOURCE_DIRECTORY(${PROJECT_SOURCE_DIR}/src/atlasstyle ATLASSRC)
AUX_SOURCE_DIRECTORY(${PROJECT_SOURCE_DIR}/src/external EXTERNALSRC)
AUX_SOURCE_DIRECTORY(${PROJECT_SOURCE_DIR}/src/fcnc FCNCSRC)
add_library(AtlasStyle SHARED ${ATLASSRC})
add_library(External SHARED ${EXTERNALSRC})
add_library(PlotTool SHARED ${FCNCSRC})
target_link_libraries(PlotTool External AtlasStyle ${ROOT_LIBRARIES})
target_link_libraries(AtlasStyle ${ROOT_LIBRARIES})
add_executable(test_run ${PROJECT_SOURCE_DIR}/util/test.cc)
target_link_libraries(test_run External)