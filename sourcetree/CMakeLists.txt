cmake_minimum_required(VERSION 2.8)

project(pdebc C CXX)


OPTION(BUILD_DOCUMENTATION "Use Doxygen to create the HTML based API documentation" OFF)
if(BUILD_DOCUMENTATION)
	FIND_PACKAGE(Doxygen)
	if (NOT DOXYGEN_FOUND)
		MESSAGE(FATAL_ERROR "Doxygen is needed to build the documentation but I failed to find it T_T")
	endif()
	configure_file(${CMAKE_SOURCE_DIR}/../doc/Doxyfile ${CMAKE_BINARY_DIR}/doc/Doxyfile @ONLY)
	add_custom_target(doc ALL COMMAND ${DOXYGEN_EXECUTABLE} Doxyfile WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/doc SOURCES ${CMAKE_BINARY_DIR}/doc/Doxyfile)
endif()

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

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
	# using Clang
	SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -std=c++11")
	SET(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -pipe -fomit-frame-pointer -std=c++11")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
	# using GCC
	SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -std=c++11")
	SET(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -pipe -fomit-frame-pointer -std=c++11")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
	# using Intel C++
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
	# using Visual Studio C++
endif()


add_library(pdebc SHARED ${SRCS})
target_link_libraries(pdebc ${CMAKE_THREAD_LIBS_INIT})

install (TARGETS pdebc DESTINATION lib/pdebc)
install (FILES ${HEADERS} DESTINATION include/pdebc)
