file(GLOB HEADERS 
        ${CMAKE_CURRENT_SOURCE_DIR}/*.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/glad/include/glad/glad.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/gui/*.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/shader/*.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/graphics/*.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/utils/*.hpp
        )
file(GLOB SOURCES 
        ${CMAKE_CURRENT_SOURCE_DIR}/glad/*.c
        ${CMAKE_CURRENT_SOURCE_DIR}/gui/*.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/shader/*.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/graphics/*.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/utils/*.cpp
        )

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
add_library(Garuda STATIC ${SOURCES} ${HEADERS})
add_dependencies(Garuda glfw)
add_dependencies(Garuda tinyobjloader)
add_dependencies(Garuda Ramuh)
target_link_libraries(Garuda PUBLIC ${GLFW_LIBRARIES} glad Ramuh tinyobjloader glfw)
set_target_properties(Garuda PROPERTIES
        OUTPUT_NAME "Garuda"
        FOLDER "GARUDA")
target_include_directories(Garuda PUBLIC ${GARUDA_SOURCE_DIR}/garuda )

set(GARUDA_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}")
# add_library(glad STATIC "${GLAD_DIR}/glad.c")
# SET_TARGET_PROPERTIES(glad PROPERTIES LINKER_LANGUAGE CXX)
# target_include_directories(glad PRIVATE "${GLAD_DIR}/include/glad.h")