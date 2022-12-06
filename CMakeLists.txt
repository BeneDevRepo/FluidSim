cmake_minimum_required(VERSION 3.0.0)

set(CMAKE_CXX_STANDARD 20)

project(Client)
add_compile_options(-O2)
add_compile_definitions(UNICODE)

# -- Main Executable --
file(GLOB_RECURSE SOURCES "src/*.cpp") # All c++ files
file(GLOB_RECURSE MODULES_SOURCES "modules/*.cpp") # All c++ files
add_executable(Client main.cpp ${SOURCES} ${MODULES_SOURCES})

# -- Headers --
target_include_directories(Client PUBLIC "modules")
target_include_directories(Client PUBLIC "src")

# -- Libraries --
target_link_libraries(Client "ws2_32.lib")
target_link_libraries(Client "winmm.lib")
# target_link_libraries(Client "libopengl32.lib")