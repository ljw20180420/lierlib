cmake_minimum_required(VERSION 3.22.1)
project(rearrangement)

add_executable(rearrangement main.cpp)
target_compile_options(rearrangement PRIVATE -std=c++20)

install(TARGETS rearrangement RUNTIME DESTINATION bin)
