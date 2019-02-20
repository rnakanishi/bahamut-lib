file(GLOB HEADERS 
        ${CMAKE_CURRENT_SOURCE_DIR}/*.h
        ${CMAKE_CURRENT_SOURCE_DIR}/glad/include/glad/glad.h
        ${CMAKE_CURRENT_SOURCE_DIR}/gui/*.h
        )
file(GLOB SOURCES 
        ${CMAKE_CURRENT_SOURCE_DIR}/glad/*.c
        ${CMAKE_CURRENT_SOURCE_DIR}/gui/*.cpp
        )

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
add_library(Garuda STATIC ${SOURCES} ${HEADERS})
target_link_libraries(Garuda PUBLIC ${GLFW_LIBRARIES} glad)
set_target_properties(Garuda PROPERTIES
        OUTPUT_NAME "Garuda"
        FOLDER "GARUDA")
target_include_directories(Garuda PUBLIC ${GARUDA_SOURCE_DIR}/garuda)

set(GARUDA_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}")
# add_library(glad STATIC "${GLAD_DIR}/glad.c")
# SET_TARGET_PROPERTIES(glad PROPERTIES LINKER_LANGUAGE CXX)
# target_include_directories(glad PRIVATE "${GLAD_DIR}/include/glad.h")