add_executable (semi_lagrangian_adv_2d_case.exe "${PROJECT_SOURCE_DIR}/semi_lagrangian/semi_lagrangian_adv_2d_case.f90")
target_link_libraries (semi_lagrangian_adv_2d_case.exe netcdf netcdff)