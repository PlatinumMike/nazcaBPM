cmake_minimum_required(VERSION 3.28)
project(testSuite)

set(CMAKE_CXX_STANDARD 20)

#todo: merge this with the CMakeLists one level up so that you can re-use the precompiled source files from the main project. But, make sure that building the tests is an opt-in, because the average user won't touch this.
add_executable(testSuite.x test.cpp ../src/TriDiag.cpp ../src/TriDiag.h)