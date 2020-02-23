program linear_swm_a_grid

  use linear_swm_common_mod

  implicit none

  call init_swm()

  ! Set initial conditions.
  call set_ic_h_spike(old)
  call output(old)

  dt = dx / sqrt(g * h0) * 0.5_rp
  do while (time_step < 20)
    call run_swm(old, new)
    call output(new)
    call time_advance()
  end do

contains

  subroutine run_swm(old, new)

    integer, intent(in) :: old, new

    integer i, j

    !  _____________________________
    ! |         |         |         |
    ! |         |         |         |
    ! |   huv   |   huv   |   huv   |  j+1
    ! |         |         |         |
    ! |_________|_________|_________|
    ! |         |         |         |
    ! |         |         |         |
    ! |   huv   |   huv   |   huv   |   j
    ! |         |         |         |
    ! |_________|_________|_________|
    ! |         |         |         |
    ! |         |         |         |
    ! |   huv   |   huv   |   huv   |  j-1
    ! |         |         |         |
    ! |_________|_________|_________|
    !
    !     i-1        i        i+1               <- full x grids
    !                                   ^
    !                                   |
    !                          full y grids

    do j = 1, ny
      do i = 1, nx
        du(i,j,old) = &
          ! Coriolis force
          f * v(i,j,old) &
          ! Pressure gradient force
         -g * (h(i+1,j,old) - h(i-1,j,old)) / (2 * dx)
        dv(i,j,old) = &
          ! Coriolis force
         -f * u(i,j,old) &
          ! Pressure gradient force
         -g * (h(i,j+1,old) - h(i,j-1,old)) / (2 * dy)
        dh(i,j,old) = &
          ! Mass flux gradient along x direction
         -h0 * (u(i+1,j,old) - u(i-1,j,old)) / (2 * dx) &
          ! Mass flux gradient along y direction
         -h0 * (v(i,j+1,old) - v(i,j-1,old)) / (2 * dy)
      end do
    end do

    do j = 1, ny
      do i = 1, nx
        u(i,j,new) = u(i,j,old) + dt * du(i,j,old)
        v(i,j,new) = v(i,j,old) + dt * dv(i,j,old)
        h(i,j,new) = h(i,j,old) + dt * dh(i,j,old)
      end do
    end do
    call apply_bc(new)

  end subroutine run_swm

  subroutine output(k)

    use netcdf

    integer, intent(in) :: k

    integer ncid, ierr
    integer x_dimid, y_dimid, time_dimid
    integer xi_dimid, yi_dimid
    integer x_varid, y_varid, xi_varid, yi_varid, time_varid
    integer u_varid, v_varid, h_varid

    character(256) file_name

    write(file_name, "('linear_swm_a_grid.', I3.3, '.nc')") time_step

    ierr = nf90_create(file_name, nf90_clobber, ncid)

    ierr = nf90_def_dim(ncid, 'time', nf90_unlimited, time_dimid)

    ierr = nf90_def_var(ncid, 'time', nf90_int, [time_dimid], time_varid)

    ierr = nf90_def_dim(ncid, 'x', nx, x_dimid)

    ierr = nf90_def_var(ncid, 'x', nf90_double, [x_dimid], x_varid)

    ierr = nf90_def_dim(ncid, 'xi', nx, xi_dimid)

    ierr = nf90_def_var(ncid, 'xi', nf90_double, [xi_dimid], xi_varid)

    ierr = nf90_def_dim(ncid, 'y', nx, y_dimid)

    ierr = nf90_def_var(ncid, 'y', nf90_double, [y_dimid], y_varid)

    ierr = nf90_def_dim(ncid, 'yi', ny, yi_dimid)

    ierr = nf90_def_var(ncid, 'yi', nf90_double, [yi_dimid], yi_varid)

    ierr = nf90_def_var(ncid, 'u', nf90_double, [x_dimid,y_dimid,time_dimid], u_varid)

    ierr = nf90_def_var(ncid, 'v', nf90_double, [x_dimid,y_dimid,time_dimid], v_varid)

    ierr = nf90_def_var(ncid, 'h', nf90_double, [x_dimid,y_dimid,time_dimid], h_varid)

    ierr = nf90_enddef(ncid)

    ierr = nf90_put_var(ncid, time_varid, time_step)

    ierr = nf90_put_var(ncid, x_varid, full_x(1:nx))

    ierr = nf90_put_var(ncid, xi_varid, half_x(1:nx))

    ierr = nf90_put_var(ncid, y_varid, full_y(1:ny))

    ierr = nf90_put_var(ncid, yi_varid, half_y(1:ny))

    ierr = nf90_put_var(ncid, u_varid, u(1:nx,1:ny,k))

    ierr = nf90_put_var(ncid, v_varid, v(1:nx,1:ny,k))

    ierr = nf90_put_var(ncid, h_varid, h(1:nx,1:ny,k))

    ierr = nf90_close(ncid)

  end subroutine output

end program linear_swm_a_grid
