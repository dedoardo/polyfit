cmake_minimum_required(VERSION 3.0)

project(polyvec-apps)

message("----------------------------------------------------------------------")
message("build applications")

# Todo this should be part of polyvec-lib (or not)
file(GLOB_RECURSE POLYVEC_DEV_FILES ${CMAKE_SOURCE_DIR}/source/dev/*.cpp)
message("polyvec files")
foreach(FILE_SRC ${POLYVEC_DEV_FILES})
	message("\t${FILE_SRC}")
endforeach(FILE_SRC)

# Enviromental variables set by the parent
message("POLYVEC_LIB_INCLUDES")
foreach(DIR ${POLYVEC_LIB_INCLUDES})
	message("\t${DIR}")
endforeach()

message("POLYVEC_LIB_LIBRARIES")
foreach(LIB ${POLYVEC_LIB_LIBRARIES})
	message("\t${LIB}")
endforeach()


# Adding all applications in dev/
# file(GLOB_RECURSE DEV_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/dev/*.cpp)
# Adding only the main application
set(DEV_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/polyfit.cpp)

# ==========================================================
# Build each file in apps/dev into its own executable
# ==========================================================
if (1)
	foreach(FILE_SRC ${DEV_SOURCES})
		get_filename_component(FILE_EXE ${FILE_SRC} NAME_WE)
		message("build file: ${FILE_SRC} -> ${FILE_EXE}")

		add_executable(${FILE_EXE} ${FILE_SRC} ${POLYVEC_DEV_FILES})
		target_include_directories(${FILE_EXE} SYSTEM PUBLIC
            "${CMAKE_SOURCE_DIR}/include"
			"${CMAKE_SOURCE_DIR}/dependencies"
            "${CMAKE_SOURCE_DIR}/dependencies/eigen"
            "${CMAKE_SOURCE_DIR}/dependencies/cairo"
            "${CMAKE_SOURCE_DIR}/dependencies/stb"
		)

		if (WIN32)
		target_include_directories(${FILE_EXE} SYSTEM PUBLIC
            "${CMAKE_SOURCE_DIR}/dependencies/dirent"
		)
		endif()

		target_link_libraries(${FILE_EXE} PUBLIC polyvec-lib)
		if (UNIX)
			target_link_libraries(${FILE_EXE} PUBLIC cairo)
		else()
			target_link_libraries(${FILE_EXE} PUBLIC "${CMAKE_SOURCE_DIR}/dependencies/cairo/x64/cairo.lib")
		endif()


		set_target_properties(${FILE_EXE}
			PROPERTIES
			RUNTIME_OUTPUT_DIRECTORY ${POLYVEC_BIN_DIR})

		target_compile_definitions(${FILE_EXE} PUBLIC ${POLYVEC_DEFINITIONS})
	endforeach(FILE_SRC)
endif()

# ==========================================================
# All the files in the subdirectory are linked together with tests.cpp
# tests.cpp contains the application entry points which imports and executes
# the tests.
# ==========================================================
if (0)
	file(GLOB_RECURSE TEST_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/tests/*.cpp)

	message("----------------------------------------------------------------------")
	message("build tests")

	foreach(FILE_SRC ${TEST_SOURCES})
		message("\t${FILE_SRC}")
	endforeach(FILE_SRC)

	add_executable(tests ${TEST_SOURCES} pvec-test.cpp ${POLYVEC_DEV_FILES})
	target_include_directories(tests SYSTEM PUBLIC ${POLYVEC_LIB_INCLUDES})
	target_link_libraries(tests ${POLYVEC_LIB} ${POLYVEC_EXTERNAL_LIB})
	target_compile_definitions(tests PUBLIC
		"-DCYCPPUTIL_TEST_MODULES_PATH=\"${CMAKE_SOURCE_DIR}/tests/cycpputil/\""
		"-DPOLYVEC_TEST_PATH=\"${CMAKE_SOURCE_DIR}/tests/\""
	)

	set_target_properties(tests
		PROPERTIES
		RUNTIME_OUTPUT_DIRECTORY ${POLYVEC_BIN_DIR}
		VS_DEBUGGER_WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
	)

	target_compile_definitions(${FILE_EXE} PUBLIC ${POLYVEC_DEFINITIONS})

endif()
