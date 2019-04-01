
add_executable(levelset main.cpp levelset_fluid.cpp)
target_link_libraries(levelset Ramuh)

add_executable(levelset3 main3.cpp levelset_fluid.cpp)
target_link_libraries(levelset3 Ramuh)
