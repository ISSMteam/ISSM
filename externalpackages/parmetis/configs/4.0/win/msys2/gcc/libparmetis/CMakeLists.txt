# Include directories for library code.
include_directories(.)
# Find sources.
file(GLOB parmetis_sources *.c)
# Create libparmetis
add_library(parmetis ${ParMETIS_LIBRARY_TYPE} ${parmetis_sources})

# Link with metis and MPI libraries.
set(CMAKE_FIND_LIBRARY_PREFIXES "lib" "msys-" "")
set(CMAKE_FIND_LIBRARY_SUFFIXES ".dll" ".so")
find_library(METIS metis PATHS ${METIS_PATH}/lib)
if(METIS)
	message("-- Found user-supplied METIS: ${METIS}")
	target_link_libraries(parmetis ${METIS})
else(METIS)
	target_link_libraries(parmetis metis)
endif()

if(LOCAL_MSMPI)
	# Add local MS-MPI libraries
	find_library(MSMPI msmpi PATHS ${MSMPI_ROOT}/lib)
	message(STATUS "MSMPI: ${MSMPI}")
	target_link_libraries(parmetis ${MSMPI})
else(LOCAL_MSMPI)
	target_link_libraries(parmetis ${MPI_LIBRARIES})
	set_target_properties(parmetis PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()

install(TARGETS parmetis
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION lib
  ARCHIVE DESTINATION lib)
