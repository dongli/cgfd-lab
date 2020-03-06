module adv_1d_common_mod

  use netcdf

  implicit none

  integer :: time_step = 0, old = 1, new = 2

contains

  subroutine advance_time()

    integer tmp

    ! Change time indices.
    tmp = old; old = new; new = tmp
    time_step = time_step + 1

  end subroutine advance_time

  subroutine apply_bc(ns, nx, x)

    integer, intent(in) :: ns
    integer, intent(in) :: nx
    real, intent(inout) :: x(1-ns:nx+ns)

    integer i

    do i = 1, ns
      x(1-i) = x(nx-i+1)
      x(nx+i) = x(1+i-1)
    end do

  end subroutine apply_bc

  subroutine output(scheme, time_step, ns, nx, x, rho)

    character(*), intent(in) :: scheme
    integer, intent(in) :: time_step
    integer, intent(in) :: ns
    integer, intent(in) :: nx
    real, intent(in) :: x(nx)
    real, intent(in) :: rho(1-ns+1:nx+ns)

    character(30) file_name
    integer ncid, time_dimid, time_varid, x_dimid, x_varid, rho_varid, ierr

    write(file_name, "(A, '.', I3.3, '.', I4.4, '.nc')") scheme, nx, time_step

    ierr = NF90_CREATE(file_name, NF90_CLOBBER, ncid)
    ierr = NF90_PUT_ATT(ncid, NF90_GLOBAL, 'scheme', scheme)

    ierr = NF90_DEF_DIM(ncid, 'time', NF90_UNLIMITED, time_dimid)
    ierr = NF90_DEF_VAR(ncid, 'time', NF90_INT, [time_dimid], time_varid)
    ierr = NF90_DEF_DIM(ncid, 'x', nx, x_dimid)
    ierr = NF90_DEF_VAR(ncid, 'x', NF90_FLOAT, [x_dimid], x_varid)
    ierr = NF90_DEF_VAR(ncid, 'rho', NF90_FLOAT, [x_dimid, time_dimid], rho_varid)

    ierr = NF90_ENDDEF(ncid)

    ierr = NF90_PUT_VAR(ncid, time_varid, time_step)
    ierr = NF90_PUT_VAR(ncid, x_varid, x)
    ierr = NF90_PUT_VAR(ncid, rho_varid, rho(1:nx))

    ierr = NF90_CLOSE(ncid)

  end subroutine output

end module adv_1d_common_mod