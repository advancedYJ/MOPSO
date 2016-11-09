cmake_minimum_required(VERSION 3.5)
project(version_2)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -lpthread")

set(SOURCE_FILES
    CMakeLists.txt
    coordinaryTransform.cpp
    fileDisposal.cpp
    initialization.cpp
    main.cpp
    matrixFunction.cpp
    MOPSO.h
    MOPSOAidFunction.cpp
    MOPSOFunction.cpp
    parameterSettings.cpp)

add_executable(version_2 ${SOURCE_FILES})