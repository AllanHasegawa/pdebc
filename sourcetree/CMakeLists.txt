
project(pdebc)

set(SRCS
	main.cpp
	external/tinyxml2/tinyxml2.cpp
	SVGParser.cpp
)

SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -std=c++11")

add_executable(pdebc ${SRCS})
