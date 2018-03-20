add_executable (upwind_adv_1d_case.exe "${PROJECT_SOURCE_DIR}/upwind/upwind_adv_1d_case.f90")
target_link_libraries (upwind_adv_1d_case.exe netcdf netcdff)