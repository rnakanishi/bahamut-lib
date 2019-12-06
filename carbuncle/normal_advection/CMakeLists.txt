
file(GLOB HEADER_FILES ./*.hpp)

file( GLOB SOURCES ./*.cpp)

add_executable(normalAdvection ${SOURCES} ${HEADER_FILES})
target_include_directories(normalAdvection PRIVATE .)
target_link_libraries(normalAdvection Ramuh Leviathan pugixml)