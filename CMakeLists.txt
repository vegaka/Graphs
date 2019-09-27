cmake_minimum_required(VERSION 3.10)
project(graphs)

set(CMAKE_CXX_STANDARD 11)

add_executable(graphs)

target_include_directories(graphs PRIVATE include)

add_subdirectory(src)