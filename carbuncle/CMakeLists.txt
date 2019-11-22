
cmake_minimum_required(VERSION 3.1)
project(CARBUNCLE LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)


# # include_directories(${CMAKE_CURRENT_SOURCE_DIR})
# add_library(Carbuncle STATIC ${SOURCES} ${HEADERS})
# add_dependencies(Carbuncle Ramuh)
# add_dependencies(Carbuncle Leviathan)
# add_dependencies(Carbuncle Garuda)
# add_dependencies(Carbuncle glm)
# target_link_libraries(Carbuncle PUBLIC Ramuh Leviathan Garuda)
# set_target_properties(Carbuncle PROPERTIES
#         OUTPUT_NAME "Carbuncle"
#         FOLDER "CARBUNCLE")
# target_include_directories(Carbuncle PUBLIC ${CARBUNCLE_SOURCE_DIR} ${GLM_INCLUDE_DIRS})


add_subdirectory(dualcubes)
add_subdirectory(weno)
add_subdirectory(levelset)
add_subdirectory(particles)
add_subdirectory(enright)
add_subdirectory(heat_kernel)