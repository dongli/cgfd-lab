add_executable(ffsl_adv_1d_case.exe ffsl_adv_1d_case.f90)
target_link_libraries(ffsl_adv_1d_case.exe netcdff)

add_executable(ffsl_adv_2d_case.exe ffsl_adv_2d_case.f90)
target_link_libraries(ffsl_adv_2d_case.exe netcdff)
