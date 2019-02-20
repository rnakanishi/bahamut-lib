
message("Building ui executable")
message("glad dir: " ${GARUDA_INCLUDE_DIR})

add_executable(ui hello.cpp )
target_link_libraries(ui PUBLIC Garuda ${GLFW_LIBRARIES} glad)
target_include_directories( ui PUBLIC ${GLFW_INCLUDE_DIR} ${GLAD_INCLUDES})
# add_dependencies( ui glfw )

