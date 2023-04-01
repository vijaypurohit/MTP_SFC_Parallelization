cmake_minimum_required(VERSION 3.24)
project(SFC_Parallelization)

set(CMAKE_CXX_STANDARD 20)

include_directories(.)

add_executable(SFC_Parallelization
        main.cpp
        PhysicalGraph.h VirtualMachines.h ServiceFunctionChain.h FileReadingFunctions.h TimeCalculationFunctions.h AssignmentFunctions.h Algorithms.h)
