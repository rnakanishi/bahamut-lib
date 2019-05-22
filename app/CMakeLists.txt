
add_subdirectory(ui-test)
add_subdirectory(app_levelset)

add_executable(linearInterpolation linear.cpp)
target_link_libraries(linearInterpolation Ramuh)

add_executable(cubicInterpolation cubic.cpp)
target_link_libraries(cubicInterpolation Ramuh)

add_executable(objloader loader_test.cpp)
target_link_libraries(objloader PUBLIC Garuda ${GLFW_LIBRARIES} glad)
target_include_directories( objloader PUBLIC ${GLFW_INCLUDE_DIR} ${GLAD_INCLUDES})
