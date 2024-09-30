cmake_minimum_required(VERSION 3.28)
project(nazcaBPM)

set(CMAKE_CXX_STANDARD 20)

add_executable("${PROJECT_NAME}.x" src/main.cpp
        src/Engine.cpp
        src/Engine.h
        src/AuxiliaryFunctions.cpp
        src/AuxiliaryFunctions.h
        src/Writers.cpp
        src/Writers.h)