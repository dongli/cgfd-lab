find_package (GSL)

if (${GSL_FOUND})
	include_directories (${GSL_INCLUDE_DIRS})
	link_directories (${GSL_LIBRARY_DIRS})
	enable_language (C)
	add_executable (crank_nicolson_adv_1d_case.exe "${PROJECT_SOURCE_DIR}/crank_nicolson/crank_nicolson_adv_1d_case.f90")
	target_link_libraries (crank_nicolson_adv_1d_case.exe netcdf netcdff gsl)
endif ()