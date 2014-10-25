cmake_minimum_required(VERSION 2.8)

project(ValaisSun)

find_package(Boost 1.54.0 REQUIRED COMPONENTS
	thread
	filesystem
	system
	iostreams
)

include_directories(
	/opt/local/include
	${PROJECT_SOURCE_DIR}/gnuplot-iostream
)

add_executable(ValaisSun
	sources/main.cpp
	sources/vector3d.hpp
	sources/position.hpp
)

find_library(OPENGL_LIB
	OpenGL
)

target_link_libraries(ValaisSun
	${OPENGL_LIB}
	${Boost_LIBRARIES}
)	
