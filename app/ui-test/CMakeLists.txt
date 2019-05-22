
message("Building ui executable")

add_executable(ui hello.cpp)
target_link_libraries(ui PUBLIC Garuda ${GLFW_LIBRARIES} glad)
target_include_directories( ui PUBLIC ${GLFW_INCLUDE_DIR} ${GLAD_INCLUDES})

add_executable(texture texture.cpp)
target_link_libraries(texture PUBLIC Garuda ${GLFW_LIBRARIES} glad)
target_include_directories( texture PUBLIC ${GLFW_INCLUDE_DIR} ${GLAD_INCLUDES})

add_executable(surface surface.cpp)
target_link_libraries(surface PUBLIC Garuda ${GLFW_LIBRARIES} glad)
target_include_directories( surface PUBLIC ${GLFW_INCLUDE_DIR} ${GLAD_INCLUDES})

add_executable(camera camera_test.cpp)
target_link_libraries(camera PUBLIC Garuda ${GLFW_LIBRARIES} glad)
target_include_directories( camera PUBLIC ${GLFW_INCLUDE_DIR} ${GLAD_INCLUDES})


