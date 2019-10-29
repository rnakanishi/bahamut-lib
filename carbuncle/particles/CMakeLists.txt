
add_executable(psystem2 psystem2.cpp)
target_link_libraries(psystem2 Ramuh)

add_executable(plevelset2 plevelset2.cpp)
target_link_libraries(plevelset2 Ramuh Leviathan)

add_executable(deformation deformation.cpp)
target_link_libraries(deformation Ramuh Leviathan)

add_executable(rigid3d rigid3d.cpp)
target_link_libraries(rigid3d Ramuh Leviathan)

add_executable(deform3d deform3d.cpp)
target_link_libraries(deform3d Ramuh Leviathan)

add_executable(runConfig run_from_config.cpp ../dualcubes/dual_cubes3.cpp)
target_link_libraries(runConfig Ramuh Leviathan pugixml)
