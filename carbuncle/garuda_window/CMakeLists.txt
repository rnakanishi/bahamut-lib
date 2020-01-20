

# add_executable(garudaWindow garuda_window.cpp)
# target_link_libraries(garudaWindow Ramuh Leviathan Garuda nanogui ${NANOGUI_EXTRA_LIBS})

add_executable(nanoguiTest nanogui_test.cpp)
target_link_libraries(nanoguiTest nanogui ${NANOGUI_EXTRA_LIBS})
