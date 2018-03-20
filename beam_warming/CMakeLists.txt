add_executable (beam_warming_adv_1d_case.exe "${PROJECT_SOURCE_DIR}/beam_warming/beam_warming_adv_1d_case.f90")
target_link_libraries (beam_warming_adv_1d_case.exe netcdf netcdff)