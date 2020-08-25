cmake_minimum_required(VERSION 2.6)
project(polyvec)

add_compile_options($<$<CXX_COMPILER_ID:MSVC>:/MP>)
add_compile_options($<$<CXX_COMPILER_ID:MSVC>:/EHsc>)

# C++11
set(CMAKE_CXX_STANDARD 11)

# Multiprocessor compilation
if (WIN32)
    set(CMAKE_CXX_FLAGS /MP)
    set(CMAKE_C_FLAGS /MP)
endif()

# Making sure that find_package() works on Windows
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# name of the header and source file directories where files should be looked for 
# under root/include and root/source
set(POLYVEC_MODULES   
				"core" 
				"curve-tracer"
				"mc"
				"polygon-tracer"
				"regularity"
				"geometry"
				"image-segment"
				"visualize"
				"shortest-path"
				"utils"
				"io"
)

foreach(MODULE ${POLYVEC_MODULES})
	message("found module ${MODULE}")
endforeach(MODULE)

# Projects to be build
set(POLYVEC_BIN_DIR "${CMAKE_CURRENT_SOURCE_DIR}" CACHE STRING "Where to write the application files.")
set(POLYVEC_LIB_DIR "${CMAKE_CURRENT_SOURCE_DIR}" CACHE STRING "Where to write the library files.")

# Only 64-bit are supported
if(${CMAKE_SIZEOF_VOID_P} MATCHES "4")
	message(FATAL_ERROR "The library might work, but alignment on Eigen containers is not guaranteed")
endif()

# Globbing all the include and source files for each module
message("finding header and source files")
file(GLOB POLYVEC_SOURCES include/polyvec/*.hpp source/polyvec/*.cpp)
message("old sources: ${POLYVEC_SOURCES}")
foreach(MODULE ${POLYVEC_MODULES})
	set(MODULE_INCLUDES_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include/polyvec/${MODULE})
	set(MODULE_SOURCES_DIR ${CMAKE_CURRENT_SOURCE_DIR}/source/polyvec/${MODULE})

	message("${MODULE} headers in ${MODULE_INCLUDES_DIR}")
	FILE(GLOB_RECURSE MODULE_INCLUDES ${MODULE_INCLUDES_DIR}/*.hpp)
	SET(POLYVEC_SOURCES ${POLYVEC_SOURCES} ${MODULE_INCLUDES})
	foreach(INC ${MODULE_INCLUDES})
		message("\t${INC}")
	endforeach(INC)

	message("${MODULE} sources in ${MODULE_SOURCES_DIR}")
	FILE(GLOB_RECURSE MODULE_SOURCES ${MODULE_SOURCES_DIR}/*.cpp)
	SET(POLYVEC_SOURCES ${POLYVEC_SOURCES} ${MODULE_SOURCES})
	foreach(SRC ${MODULE_SOURCES})
		message("\t${SRC}")
	endforeach(SRC)
endforeach(MODULE)

# Main library will be written in lib/
add_library(polyvec-lib ${POLYVEC_SOURCES})
set_target_properties(polyvec-lib
	PROPERTIES
	ARCHIVE_OUTPUT_DIRECTORY ${POLYVEC_LIB_DIR})

# Include directories for the vectorization library
target_include_directories(polyvec-lib PUBLIC 
    "${CMAKE_CURRENT_SOURCE_DIR}/include" 
)

# In order to link against cairo on Windows using MSVC, a x64 .lib/.dll binaries are bundled in the repository
# On Linux you can install it with 'apt install libcairo2-dev'.
# Apparently cairo doesn't provide the cmake configuration files: https://github.com/Microsoft/vcpkg/issues/2563 
if (${CMAKE_CXX_COMPILER_ID} STREQUAL "MSVC")
	target_compile_definitions(polyvec-lib PUBLIC "-DPOLYVEC_USE_CAIRO=1")
	target_link_libraries(polyvec-lib ${CMAKE_CURRENT_SOURCE_DIR}/dependencies/cairo/x64/cairo.lib)
	target_include_directories(polyvec-lib SYSTEM PUBLIC 
        "${CMAKE_CURRENT_SOURCE_DIR}/dependencies/cairo"
    )
	file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/dependencies/cairo/x64/cairo.dll DESTINATION ${CMAKE_CURRENT_SOURCE_DIR})
elseif((${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU") OR
		(${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang"))

	find_library(
		Cairo_LIBRARIES 
			cairo
		HINTS
			"/usr/lib"
			"/usr/lib/x86_64-linux-gnu"
			"/usr/local/lib"
	)

	find_path(
		Cairo_INCLUDE_DIRS 
			cairo/cairo.h
		HINTS
			"/usr/include"
			"/usr/local/include"
	)

	target_include_directories(polyvec-lib SYSTEM PUBLIC ${Cairo_INCLUDE_DIRS})
	target_link_libraries(polyvec-lib PUBLIC ${Cairo_LIBRARIES})
else()
	message(FATAL_ERROR "Unrecognized compiler ${CMAKE_CXX_COMPILER_ID}")
endif()

# Random forest 
# Great implementation
# https://github.com/bjoern-andres/random-forest
message("dependency: Andres' Random Forest")
target_include_directories(polyvec-lib SYSTEM PUBLIC 
	"${CMAKE_CURRENT_SOURCE_DIR}/dependencies/random-forest/include"
)

# Eigen3.3 +
# Todo: provide a FindEigen3.cmake
target_include_directories(polyvec-lib SYSTEM PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/dependencies/eigen")

if (WIN32)
	target_compile_definitions(polyvec-lib PUBLIC "-D_ENABLE_EXTENDED_ALIGNED_STORAGE")
endif()

if (MSVC)
	target_compile_definitions(polyvec-lib PUBLIC "-D_CRT_SECURE_NO_WARNINGS")
endif()

# more header-only dependencies
target_include_directories(polyvec-lib SYSTEM PUBLIC
    # Used to list files in directory
    "${CMAKE_SOURCE_DIR}/dependencies/dirent"

    # Used to read raster images
    "${CMAKE_SOURCE_DIR}/dependencies/stb"
)


# Building executables depending on it
add_subdirectory(apps)
