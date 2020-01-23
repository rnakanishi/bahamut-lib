
add_executable(dualCubes3 main.cpp dual_cubes3.cpp)
target_link_libraries(dualCubes3 Ramuh Leviathan)

add_executable(dualCubes2 main2d.cpp dual_cubes2.cpp)
target_link_libraries(dualCubes2 Ramuh Leviathan)

add_executable(normalParticle2 main_normalParticles.cpp dual_cubes2.cpp normal_particles2.cpp)
target_link_libraries(normalParticle2 Ramuh Leviathan)

add_executable(tangentParticle2 main_tangentParticles.cpp dual_cubes2.cpp tangent_particles2.cpp )
target_link_libraries(tangentParticle2 Ramuh Leviathan)

add_executable(fusingParticle2 main_fuseParticles.cpp dual_cubes2.cpp tangent_particles2.cpp )
target_link_libraries(fusingParticle2 Ramuh Leviathan)

add_executable(extractSurface extract.cpp dual_cubes3.cpp)
target_link_libraries(extractSurface Ramuh Leviathan pugixml)

add_executable(tangentParticle3 main_tangentParticles3.cpp dual_cubes3.cpp tangent_particles3.cpp )
target_link_libraries(tangentParticle3 Ramuh Leviathan)