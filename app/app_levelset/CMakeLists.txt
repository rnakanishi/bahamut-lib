
add_executable(levelset main.cpp ../levelset_fluid2.cpp)
add_dependencies(levelset pugixml)
target_include_directories(levelset PUBLIC ..)
target_link_libraries(levelset Ramuh pugixml)

add_executable(levelset3 main3.cpp ../levelset_fluid3.cpp)
add_dependencies(levelset3 pugixml)
target_include_directories(levelset3 PUBLIC ..)
target_link_libraries(levelset3 Ramuh pugixml)

add_executable(triangulate triangulate.cpp ../levelset_fluid3.cpp)
target_link_libraries(triangulate Ramuh)

add_executable(marching marching.cpp ../levelset_fluid3.cpp)
target_link_libraries(marching Ramuh)
