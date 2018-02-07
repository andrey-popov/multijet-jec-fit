cmake_minimum_required(VERSION 2.8.12)
project (jecfit)


set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/lib")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/bin")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -Wall -Wextra -O2")
include_directories("${PROJECT_SOURCE_DIR}/include")


find_package(ROOT REQUIRED COMPONENTS Minuit2)
include_directories(${ROOT_INCLUDE_DIRS})

find_package(Boost 1.34.0 REQUIRED COMPONENTS program_options)
include_directories(${Boost_INCLUDE_DIR})


add_subdirectory(src)
add_subdirectory(prog)
add_subdirectory(tests)
