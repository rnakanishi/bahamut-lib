cmake_minimum_required(VERSION 3.1)
project(BAHAMUT 
        LANGUAGES CXX )
set(CMAKE_CXX_STANDARD 17)

include(ExternalProject)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/external")
include( glfw )
include( glm )

add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/external/glad" "external/glad")
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/external/pugixml" "external/pugixml")
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/external/tinyobjloader" 
                                                    "external/tinyobjloader")
set(GLAD_INCLUDES "${CMAKE_CURRENT_SOURCE_DIR}/external/glad/include")

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}\
                                 ${OpenMP_EXE_LINKER_FLAGS}")
endif()

set(EIGEN3_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/external/eigen3)
include_directories(${EIGEN3_INCLUDE_DIR})

add_subdirectory(ramuh)
add_subdirectory(garuda)
add_subdirectory(app)
