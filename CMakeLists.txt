cmake_minimum_required(VERSION 3.1)
project(BAHAMUT 
        LANGUAGES CXX CUDA)
set(CMAKE_CXX_STANDARD 17)

include(ExternalProject)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/external")
include( glfw )

add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/external/glad" "external/glad")
set(GLAD_INCLUDES "${CMAKE_CURRENT_SOURCE_DIR}/external/glad/include")

add_subdirectory(ramuh)
add_subdirectory(garuda)
add_subdirectory(app)
