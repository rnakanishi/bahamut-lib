cmake_minimum_required(VERSION 3.1)
project(BAHAMUT)
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

find_package(Eigen3 REQUIRED)
if(Eigen3_FOUND)
message("--Found Eigen3")
else()
    message("-- Building local Eigen3")
    set(EIGEN3_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/external/eigen3)
    include_directories(${EIGEN3_INCLUDE_DIR})
endif()

add_subdirectory(ramuh) # Geometric structures lib
add_subdirectory(garuda) # GUI and interface
add_subdirectory(leviathan) # Incompressible fluid simulation
add_subdirectory(carbuncle) # Incompressible fluid simulation
add_subdirectory(app)
