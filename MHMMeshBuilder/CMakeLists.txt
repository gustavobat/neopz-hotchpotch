cmake_minimum_required(VERSION 3.8...3.15)

project (NewMHMMeshBuilder LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Finding NeoPZ package
find_package(PZ REQUIRED)

add_executable(MHMMeshBuilder main.cpp)

target_link_libraries(MHMMeshBuilder pz)
target_include_directories(MHMMeshBuilder PUBLIC ${PZ_INCLUDE_DIRS})
