cmake_minimum_required(VERSION 3.8...3.15)

project (TestShape LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Finding NeoPZ package
find_package(PZ REQUIRED)

add_executable(TestShape TestShape.cpp)

target_link_libraries(TestShape pz)
target_include_directories(TestShape PUBLIC ${PZ_INCLUDE_DIRS})
