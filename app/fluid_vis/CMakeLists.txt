
add_executable(fluidVisualizer fluid_vis.cpp)
target_link_libraries(fluidVisualizer PUBLIC Garuda ${GLFW_LIBRARIES} glad)
target_include_directories( fluidVisualizer PUBLIC ${GLFW_INCLUDE_DIR} ${GLAD_INCLUDES})
