cmake_minimum_required(VERSION 2.8)

project(pdebc_bezier_fitting)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}")

find_package(LibPDEBC REQUIRED)

include_directories(${LIBPDEBC_INCLUDE_DIR})

set(SRCS
	bezier_fitting.cpp
)

SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -std=c++11 -stdlib=libc++")
SET(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -pipe -fomit-frame-pointer -std=c++11 -stdlib=libc++")
add_executable(bezier_fitting ${SRCS})
target_link_libraries(bezier_fitting ${LIBPDEBC_LIBRARY})
