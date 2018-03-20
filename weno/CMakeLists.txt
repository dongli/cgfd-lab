add_executable (weno_adv_1d_case.exe "${PROJECT_SOURCE_DIR}/weno/weno_adv_1d_case.f90")
target_link_libraries (weno_adv_1d_case.exe netcdf netcdff)