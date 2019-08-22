
add_executable(wenoVolumes wenomain.cpp)
target_link_libraries(wenoVolumes Ramuh)

add_executable(wenoDifferences wenodifferences.cpp)
target_link_libraries(wenoDifferences Ramuh)

add_executable(weno2d weno2d.cpp)
target_link_libraries(weno2d Ramuh)

add_executable(wenoAdvection2d advection2d.cpp)
target_link_libraries(wenoAdvection2d Ramuh)

add_executable(wenoAdvection1d advection1d.cpp)
target_link_libraries(wenoAdvection1d Ramuh)

add_executable(wenoRotation2d rotation2d.cpp)
target_link_libraries(wenoRotation2d Ramuh)