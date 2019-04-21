
message("Building ui executable")

add_executable(ui hello.cpp)
target_link_libraries(ui PUBLIC Garuda ${GLFW_LIBRARIES} glad)
target_include_directories( ui PUBLIC ${GLFW_INCLUDE_DIR} ${GLAD_INCLUDES})

add_executable(texture texture.cpp)
target_link_libraries(texture PUBLIC Garuda ${GLFW_LIBRARIES} glad)
target_include_directories( texture PUBLIC ${GLFW_INCLUDE_DIR} ${GLAD_INCLUDES})


