
add_executable(dualCubes3 main.cpp dual_cubes3.cpp)
target_link_libraries(dualCubes3 Ramuh Leviathan)

add_executable(dualCubes2 main2d.cpp dual_cubes2.cpp)
target_link_libraries(dualCubes2 Ramuh Leviathan)

add_executable(extractSurface extract.cpp dual_cubes3.cpp)
target_link_libraries(extractSurface Ramuh Leviathan pugixml)
