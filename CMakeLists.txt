cmake_minimum_required(VERSION 3.20)
project(nazcaBPM)

set(CMAKE_CXX_STANDARD 20)

find_package(Boost 1.70 REQUIRED)
find_package(HDF5 1.10 REQUIRED COMPONENTS CXX)

set(EXENAME "${PROJECT_NAME}.x")
add_executable(${EXENAME} src/main.cpp
        src/Engine.cpp
        src/Engine.h
        src/AuxiliaryFunctions.cpp
        src/AuxiliaryFunctions.h
        src/TriDiag.h
        src/IO/Readers.cpp
        src/IO/Readers.h
        src/Parameters.h
        src/Geometry.cpp
        src/Geometry.h
        src/Polygon.cpp
        src/Polygon.h
        src/IO/hdf_writer.h
        src/GridInterpolator.h
        src/Solver.cpp
        src/Solver.h
        src/PML.cpp
        src/PML.h
        src/OperatorSuite.cpp
        src/OperatorSuite.h
        src/RectangularGrid.cpp
        src/RectangularGrid.h
        src/IndexMonitor.cpp
        src/IndexMonitor.h
        src/FieldMonitor.cpp
        src/FieldMonitor.h
        src/Port.cpp
        src/Port.h
        src/ModeSolver.cpp
        src/ModeSolver.h
        src/BpmSolver.cpp
        src/BpmSolver.h)

target_link_libraries(${EXENAME} PRIVATE HDF5::HDF5)
target_link_libraries(${EXENAME} PRIVATE Boost::boost)

set(ENABLE_TEST 1)

if (ENABLE_TEST)
    add_executable(testSuite.x test/test.cpp src/TriDiag.h src/Polygon.cpp src/Polygon.h src/AuxiliaryFunctions.cpp src/AuxiliaryFunctions.h src/Geometry.cpp src/Geometry.h
            src/GridInterpolator.h
            src/PML.cpp
            src/PML.h
            src/OperatorSuite.cpp
            src/OperatorSuite.h
    )
endif ()

