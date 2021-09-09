add_subdirectory(io)
add_subdirectory(geometry/octree)

add_executable(nebula_cpu_vsc cpu_vsc.cpp
		$<TARGET_OBJECTS:io>)
target_link_libraries(nebula_cpu_vsc ${LIBRARIES} Threads::Threads)