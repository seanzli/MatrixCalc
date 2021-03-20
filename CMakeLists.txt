CMAKE_MINIMUM_REQUIRED(VERSION 3.10)
PROJECT(matrixText)

#include 
include_directories(${PROJECT_SOURCE_DIR}/${PROJECT_NAME})

# src
add_subdirectory(./src)
aux_source_directory(./src SRC_DIR)
add_executable(${PROJECT_NAME} ${SRC_DIR})
target_link_libraries(${PROJECT_NAME} src_lib)