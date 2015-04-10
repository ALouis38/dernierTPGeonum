# CMake entry point
cmake_minimum_required (VERSION 2.6)
project(M1_INFO)

include_directories(src/)
include_directories(
	./external/glm-0.9.4.0/
)


add_definitions(
	-DTW_STATIC
	-DTW_NO_LIB_PRAGMA
	-D_CRT_SECURE_NO_WARNINGS
)

file(
	GLOB_RECURSE
	source_files
        src/*.cpp
        src/*.h
)

add_executable(projet_m1_info ${source_files} )

