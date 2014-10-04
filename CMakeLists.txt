cmake_minimum_required(VERSION 2.8)

project(ValaisSun)

include_directories(
	/opt/local/include
)

add_executable(ValaisSun
	sources/main.cpp
)

find_library(OPENGL_LIB
	OpenGL
)

target_link_libraries(ValaisSun
	${OPENGL_LIB}
)	
