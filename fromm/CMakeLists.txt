add_executable (fromm_adv_1d_case.exe "${PROJECT_SOURCE_DIR}/fromm/fromm_adv_1d_case.f90")
target_link_libraries (fromm_adv_1d_case.exe netcdf netcdff)