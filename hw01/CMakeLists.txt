cmake_minimum_required(VERSION 3.22)
project(OSY1)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_FLAGS "-Wall -pedantic -fsanitize=thread")
link_directories(x86_64-linux-gnu)
link_libraries(pthread progtest_solver)

add_executable(OSY1 solution-megi.cpp sample_tester.cpp)
