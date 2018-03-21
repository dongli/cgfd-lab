! This is a 1D advection example using square initial condition and periodic
! boundary condition for Fromm finite difference scheme.
!
! Li Dong <dongli@lasg.iap.ac.cn>
!
! - 2018-03-20: Initial creation.

program fromm_adv_1d_case

  use netcdf

  implicit none

  real, allocatable :: x(:)         ! Cell center coordinates
  real, allocatable :: rho(:,:)     ! Tracer density being advected at cell centers
  real, allocatable :: flux(:)      ! Flux at cell interfaces
  real dx                           ! Cell interval
  real dt                           ! Time step size
  integer nx                        ! Cell number
  integer nt                        ! Integration time step number
  logical :: use_rk3 = .false.

  real :: u = 0.005                 ! Advection speed
  real coef                         ! dt / dx
  integer, parameter :: ns = 1      ! Stencil width

  integer i, time_step, old, new
  character(256) namelist_path
  logical is_exist

  namelist /params/ nx, nt, dx, dt, use_rk3, u

  call get_command_argument(1, namelist_path)
  inquire(file=namelist_path, exist=is_exist)
  if (.not. is_exist) then
    write(*, *) '[Error]: You need set the namelist path in command line!'
    stop 1
  end if

  open(10, file=namelist_path)
  read(10, nml=params)
  close(10)

  allocate(x(nx))
  allocate(rho(1-ns:nx+ns,2))
  allocate(flux(1:nx+1))

  ! Set mesh grid coordinates.
  dx = 1.0d0 / nx
  do i = 1, nx
    x(i) = (i - 1) * dx
  end do

  ! Set initial condition.
  old = 1; new = 2
  do i = 1, nx
    if (x(i) >= 0.05 .and. x(i) <= 0.3) then
      rho(i,old) = 1.0d0
    else
      rho(i,old) = 0.0d0
    end if
  end do
  call full_boundary_condition(rho(:,old))
  call output(rho(:,old))

  ! Run integration.
  coef = dt / dx
  time_step = 0
  print *, time_step, sum(rho(1:nx,old))
  do while (time_step < nt)
    ! RK 1st stage
    call fromm(rho(:,old))
    do i = 1, nx
      rho(i,new) = rho(i,old) - coef * (flux(i+1) - flux(i))
    end do
    call full_boundary_condition(rho(:,new))
    if (use_rk3) then
      ! RK 2nd stage
      call fromm(rho(:,new))
      do i = 1, nx
       rho(i,new) = (3.0d0 * rho(i,old) + rho(i,new) - coef * (flux(i+1) - flux(i))) / 4.0d0
      end do
      call full_boundary_condition(rho(:,new))
      ! RK 3st stage
      call fromm(rho(:,new))
      do i = 1, nx
       rho(i,new) = (rho(i,old) + 2.0d0 * (rho(i,new) - coef * (flux(i+1) - flux(i)))) / 3.0d0
      end do
      call full_boundary_condition(rho(:,new))
    end if
    ! Change time indices.
    i = old; old = new; new = i
    time_step = time_step + 1
    call output(rho(:,old))
    print *, time_step, sum(rho(1:nx,old))
  end do

  deallocate(x)
  deallocate(rho)
  deallocate(flux)

contains

  subroutine full_boundary_condition(x)

    real, intent(inout) :: x(1-ns:nx+ns)

    integer i

    do i = 1, ns
      x(1-i) = x(nx-i+1)
      x(nx+i) = x(1+i-1)
    end do

  end subroutine full_boundary_condition

  subroutine half_boundary_condition(x)

    real, intent(inout) :: x(1:nx+1)

    integer i

    x(1) = x(nx+1)

  end subroutine half_boundary_condition

  subroutine fromm(rho)

    real, intent(in) :: rho(1-ns:nx+ns)

    real lw, bw
    integer i

    do i = 1, nx
      lw = 0.5d0 * (u * (rho(i+1) + rho(i)) - coef * u**2 * (rho(i+1) - rho(i)))
      if (u >= 0) then
        bw = 0.5d0 * (u * (3 * rho(i  ) - rho(i-1)) - coef * u**2 * (rho(i  ) - rho(i-1)))
      else
        bw = 0.5d0 * (u * (3 * rho(i+1) - rho(i+2)) - coef * u**2 * (rho(i+2) - rho(i+1)))
      end if
      flux(i+1) = 0.5d0 * (lw + bw)
    end do
    call half_boundary_condition(flux)

  end subroutine fromm

  subroutine output(rho)

    real, intent(in) :: rho(1-ns:nx+ns)

    character(30) file_name
    integer file_id, time_dim_id, time_var_id, x_dim_id, x_var_id, rho_var_id, ierr

    write(file_name, "('fromm.', I3.3, '.nc')") time_step

    ierr = nf90_create(file_name, nf90_clobber, file_id)
    if (ierr /= nf90_noerr) then
      write(*, *) '[Error]: Failed to create output file!'
      stop 1
    end if

    ierr = nf90_def_dim(file_id, 'time', nf90_unlimited, time_dim_id)
    if (ierr /= nf90_noerr) then
      write(*, *) '[Error]: Failed to define time dimension!'
      stop 1
    end if

    ierr = nf90_def_var(file_id, 'time', nf90_int, [time_dim_id], time_var_id)
    if (ierr /= nf90_noerr) then
      write(*, *) '[Error]: Failed to define time variable!'
      stop 1
    end if

    ierr = nf90_def_dim(file_id, 'x', nx, x_dim_id)
    if (ierr /= nf90_noerr) then
      write(*, *) '[Error]: Failed to define x dimension!'
      stop 1
    end if

    ierr = nf90_def_var(file_id, 'x', nf90_float, [x_dim_id], x_var_id)
    if (ierr /= nf90_noerr) then
      write(*, *) '[Error]: Failed to define x variable!'
      stop 1
    end if

    ierr = nf90_def_var(file_id, 'rho', nf90_float, [x_dim_id, time_dim_id], rho_var_id)
    if (ierr /= nf90_noerr) then
      write(*, *) '[Error]: Failed to define rho variable!'
      stop 1
    end if

    ierr = nf90_enddef(file_id)

    ierr = nf90_put_var(file_id, time_var_id, time_step)
    if (ierr /= nf90_noerr) then
      write(*, *) '[Error]: Failed to write time variable!'
      stop 1
    end if

    ierr = nf90_put_var(file_id, x_var_id, x)
    if (ierr /= nf90_noerr) then
      write(*, *) '[Error]: Failed to write x variable!'
      stop 1
    end if

    ierr = nf90_put_var(file_id, rho_var_id, rho(1:nx))
    if (ierr /= nf90_noerr) then
      write(*, *) '[Error]: Failed to write rho variable!'
      stop 1
    end if

    ierr = nf90_close(file_id)
    if (ierr /= nf90_noerr) then
      write(*, *) '[Error]: Failed to close file!'
      stop 1
    end if

  end subroutine output

end program fromm_adv_1d_case
