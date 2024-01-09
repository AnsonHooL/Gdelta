cmake_minimum_required(VERSION 3.10)
project(Gdelta_Git)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 ")

include_directories(.)


add_executable(gdelta main.cpp gdelta.cpp)
