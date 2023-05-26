cmake_minimum_required(VERSION 3.22.1)
project(rearrangement)

add_executable(rearrangement main.cpp)
# set(CMAKE_BUILD_TYPE "Debug")
# target_compile_options(rearrangement PRIVATE -O0 -W -Wall -std=c++20)
set(CMAKE_BUILD_TYPE "Release")
target_compile_options(rearrangement PRIVATE -O3 -W -Wall -std=c++20 -static)

install(TARGETS rearrangement RUNTIME DESTINATION bin)
