cmake_minimum_required(VERSION 3.24)
project(SFC_Parallelization)

set(CMAKE_CXX_STANDARD 20)

include_directories(.)

add_executable(SFC_Parallelization
        main.cpp
        PhysicalGraph.h VirtualMachines.h ServiceFunctionChain.h FileReadingFunctions.h TimeCalculationFunctions.h AssignmentFunctions.h Algorithms.h)

#find_package(OpenMP)
#if (OPENMP_FOUND)
#    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
#    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
#    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
#endif()