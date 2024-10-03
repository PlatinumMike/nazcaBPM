cmake_minimum_required(VERSION 3.28)
project(nazcaBPM)

set(CMAKE_CXX_STANDARD 20)

add_executable("${PROJECT_NAME}.x" src/main.cpp
        src/Engine.cpp
        src/Engine.h
        src/AuxiliaryFunctions.cpp
        src/AuxiliaryFunctions.h
        src/IO/Writers.cpp
        src/IO/Writers.h
        src/TriDiag.cpp
        src/TriDiag.h
        src/IO/Readers.cpp
        src/IO/Readers.h
        src/Parameters.h
        src/Geometry.cpp
        src/Geometry.h
        src/Polygon.cpp
        src/Polygon.h)

#set(ENABLE_TEST 1)

if (ENABLE_TEST)
    add_executable(testSuite.x test/test.cpp src/TriDiag.cpp src/TriDiag.h src/Polygon.cpp src/Polygon.h src/AuxiliaryFunctions.cpp src/AuxiliaryFunctions.h src/Geometry.cpp src/Geometry.h)
endif ()