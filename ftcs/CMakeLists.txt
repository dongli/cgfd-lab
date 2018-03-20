add_executable (ftcs_adv_1d_case.exe "${PROJECT_SOURCE_DIR}/ftcs/ftcs_adv_1d_case.f90")
target_link_libraries (ftcs_adv_1d_case.exe netcdf netcdff)