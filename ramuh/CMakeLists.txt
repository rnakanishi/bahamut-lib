cmake_minimum_required(VERSION 3.1)
project(RAMUH 
        LANGUAGES CXX CUDA)
set(CMAKE_CXX_STANDARD 17)

file(GLOB HEADERS 
        ${CMAKE_CURRENT_SOURCE_DIR}/*.h
        ${CMAKE_CURRENT_SOURCE_DIR}/structures/*.h
        )
file(GLOB SOURCES 
        ${CMAKE_CURRENT_SOURCE_DIR}/structures/*.cpp
        )

# include_directories(${CMAKE_CURRENT_SOURCE_DIR})
add_library(Ramuh STATIC ${SOURCES} ${HEADERS})
set_target_properties(Ramuh PROPERTIES
        OUTPUT_NAME "Ramuh"
        FOLDER "RAMUH")
target_include_directories(Ramuh PUBLIC ${RAMUH_SOURCE_DIR})