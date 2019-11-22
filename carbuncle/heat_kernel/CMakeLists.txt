
file(GLOB HEADER_FILES ./*.hpp)

file( GLOB SOURCES ./*.cpp)

add_executable(heatKernel ${SOURCES} ${HEADER_FILES})
target_include_directories(heatKernel PRIVATE .)
target_link_libraries(heatKernel Ramuh Leviathan pugixml)