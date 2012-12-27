
project(pdebc)

set(SRCS
	main.cpp
	#BezierCurve.cpp
)

SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -std=c++11")

add_executable(pdebc ${SRCS})
