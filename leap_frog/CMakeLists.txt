add_executable (leap_frog_adv_1d_case.exe "${PROJECT_SOURCE_DIR}/leap_frog/leap_frog_adv_1d_case.f90")
target_link_libraries (leap_frog_adv_1d_case.exe netcdf netcdff)