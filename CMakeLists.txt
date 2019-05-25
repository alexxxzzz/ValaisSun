cmake_minimum_required(VERSION 2.8)

project(ValaisSun)

find_package(Boost 1.66.0 REQUIRED COMPONENTS
	thread
	filesystem
	system
	iostreams
	program_options
)

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -openmp")

include_directories(
	/opt/local/include
#	${PROJECT_SOURCE_DIR}/gnuplot-iostream
)

add_executable(CrunchGeoData
	sources/crunch.cpp
	sources/import_heightmap.hpp
	sources/swiss_to_lat_lon.hpp
	sources/vector3d.hpp
	sources/position.hpp
)

add_executable(ComputeHorizons
    sources/compute_horizons.cpp
    sources/height_map.cpp
    sources/height_map.hpp
    sources/swiss_to_lat_lon.hpp
    sources/vector3d.hpp
    sources/position.hpp
)

#add_executable(LightImage
#	sources/light.cpp
#	sources/swiss_to_lat_lon.hpp
#	sources/vector3d.hpp
#	sources/position.hpp
#)

find_library(OPENCV_CORE_LIB
    opencv_core
)
find_library(OPENCV_IMGPROC_LIB
    opencv_imgproc
)

find_library(OPENCV_HIGHGUI_LIB
    opencv_highgui
)

find_library(OPENGL_LIB
	OpenGL
)

target_link_libraries(ComputeHorizons
    ${Boost_LIBRARIES}
)

target_link_libraries(CrunchGeoData
	${OPENGL_LIB}
	${Boost_LIBRARIES}
)	

#target_link_libraries(LightImage
#	${OPENGL_LIB}
#	${Boost_LIBRARIES}
#    ${OPENCV_CORE_LIB}
#    ${OPENCV_IMGPROC_LIB}
#    ${OPENCV_HIGHGUI_LIB}
#)
