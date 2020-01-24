
file(GLOB HEADER_FILES ./*.hpp)

file( GLOB SOURCES ./*.cpp)

add_executable(enright ${SOURCES} ../dualcubes/dual_cubes3.cpp ${HEADER_FILES})
target_include_directories(enright PRIVATE .)
target_link_libraries(enright Ramuh Leviathan pugixml)