
add_executable(redistance2d redistance2.cpp)
target_link_libraries(redistance2d Ramuh Leviathan)

add_executable(redistance1d redistance1.cpp)
target_link_libraries(redistance1d Ramuh Leviathan)

add_executable(cip1d cip.cpp)
target_link_libraries(cip1d Ramuh Leviathan)
