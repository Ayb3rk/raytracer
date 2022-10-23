cmake_minimum_required(VERSION 3.21)
project(graph_hw1)

set(CMAKE_CXX_STANDARD 17)

add_executable(graph_hw1 parser.cpp parser.h ppm.cpp ppm.h raytracer.cpp tinyxml2.cpp tinyxml2.h)
