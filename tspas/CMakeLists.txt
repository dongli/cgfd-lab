add_executable (tspas_adv_1d_case.exe "${PROJECT_SOURCE_DIR}/tspas/tspas_adv_1d_case.f90")
target_link_libraries (tspas_adv_1d_case.exe netcdf netcdff)