
project(pdebc)

set(SRCS
	main.cpp
	BezierCurve.cpp
    SVGParser.cpp
    BCDESolver.cpp
	external/tinyxml2/tinyxml2.cpp
)

SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -std=c++11")

add_executable(pdebc ${SRCS})
