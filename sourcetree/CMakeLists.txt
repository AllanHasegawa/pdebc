
project(pdebc)

set(SRCS
	main.cpp
)

SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -std=c++11 -stdlib=libc++")
SET(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -pipe -fomit-frame-pointer -std=c++11 -stdlib=libc++")
add_executable(pdebc ${SRCS})
