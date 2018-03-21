add_executable (mpdata_adv_1d_case.exe "${PROJECT_SOURCE_DIR}/mpdata/mpdata_adv_1d_case.f90")
target_link_libraries (mpdata_adv_1d_case.exe netcdf netcdff)
