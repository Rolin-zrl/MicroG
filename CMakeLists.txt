cmake_minimum_required(VERSION 3.10)
project(MicroG VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_SUPPRESS_DEVELOPER_WARNINGS 1 CACHE INTERNAL "No dev warnings")

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build" FORCE)
endif()

find_package(PCL 1.12 REQUIRED)
find_package(Boost REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(OpenMP REQUIRED)
if(OPENMP_FOUND)
    message(STATUS "OpenMP found")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

set(Sophus_DIR "${CMAKE_SOURCE_DIR}/3rd/Sophus")
include_directories(${Sophus_DIR})
file(GLOB_RECURSE SOURCE_FILES "${PROJECT_SOURCE_DIR}/src/*.cpp")
file(GLOB_RECURSE HEADER_FILES "${PROJECT_SOURCE_DIR}/include/*.h" "${PROJECT_SOURCE_DIR}/include/*.hpp")


add_executable(MicroG ${SOURCE_FILES} ${HEADER_FILES} main.cpp)
target_include_directories(MicroG PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(MicroG PUBLIC
        ${PCL_LIBRARIES}
        Boost::boost
        Eigen3::Eigen
)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)