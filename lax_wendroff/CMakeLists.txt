add_executable (lax_wendroff_adv_1d_case.exe "${PROJECT_SOURCE_DIR}/lax_wendroff/lax_wendroff_adv_1d_case.f90")
target_link_libraries (lax_wendroff_adv_1d_case.exe netcdf netcdff)