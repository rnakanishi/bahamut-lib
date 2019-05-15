
add_subdirectory(ui-test)
add_subdirectory(app_levelset)

add_executable(linearInterpolation linear.cpp)
target_link_libraries(linearInterpolation Ramuh)
