cmake_minimum_required(VERSION 2.8)

project(pdebc C CXX)

FIND_PACKAGE(Threads REQUIRED)

set(SRCS
	_emptysrc.cpp
)

set(HEADERS
	SequentialDE.hpp
	BaseDE.hpp
	ThreadsDE.hpp
	ThreadsDESolver.hpp
)

SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -std=c++11")
SET(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -pipe -fomit-frame-pointer -std=c++11")
add_library(pdebc SHARED ${SRCS})
target_link_libraries(pdebc ${CMAKE_THREAD_LIBS_INIT})

install (TARGETS pdebc DESTINATION lib/pdebc)
install (FILES ${HEADERS} DESTINATION include/pdebc)
