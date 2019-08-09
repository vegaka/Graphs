cmake_minimum_required(VERSION 3.10)
project(graphs)

set(CMAKE_CXX_STANDARD 11)

include_directories(include)

add_subdirectory(src)

#add_executable(graphs ${SOURCES})
#add_executable(graphs src/main.cpp src/graph/CSR_Graph.cpp include/graph/CSR_Graph.h include/graph/Graph.h include/algorithm/SSSP.h src/algorithm/SSSP.cpp include/algorithm/SSSP.h)