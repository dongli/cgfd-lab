add_executable(weno_adv_fd_1d_case.exe
  ../adv_1d_common_mod.F90
  ../adv_1d_square_case_mod.F90
  weno_adv_fd_1d_case.f90
)
target_link_libraries(weno_adv_fd_1d_case.exe netcdff)
