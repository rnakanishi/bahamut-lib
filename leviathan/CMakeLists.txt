cmake_minimum_required(VERSION 3.1)
project(LEVIATHAN
        LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)

file(GLOB HEADERS 
        ${CMAKE_CURRENT_SOURCE_DIR}/*.h
        ${CMAKE_CURRENT_SOURCE_DIR}/fluids/*.h
        )
file(GLOB SOURCES 
        ${CMAKE_CURRENT_SOURCE_DIR}/fluids/*.cpp
        )

# include_directories(${CMAKE_CURRENT_SOURCE_DIR})
add_library(Leviathan STATIC ${SOURCES} ${HEADERS})
add_dependencies(Leviathan Ramuh)
add_dependencies(Leviathan glm)
set_target_properties(Leviathan PROPERTIES
        OUTPUT_NAME "Leviathan"
        FOLDER "LEVIATHAN")
target_include_directories(Leviathan PUBLIC ${LEVIATHAN_SOURCE_DIR} ${GLM_INCLUDE_DIRS})
