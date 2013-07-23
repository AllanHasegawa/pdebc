cmake_minimum_required(VERSION 2.8)

project(pdebc)

set(SRCS
	_emptysrc.cpp
)

set(HEADERS
	SequentialDE.hpp
	BaseDE.hpp
)

SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -std=c++11 -stdlib=libc++")
SET(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -pipe -fomit-frame-pointer -std=c++11 -stdlib=libc++")
add_library(pdebc SHARED ${SRCS})

install (TARGETS pdebc DESTINATION lib/pdebc)
install (FILES ${HEADERS} DESTINATION include/pdebc)