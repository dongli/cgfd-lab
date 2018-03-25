! This is a 1D advection example using square initial condition and periodic
! boundary condition for WENO finite difference scheme.
!
! Stencils:
!
!   For positive flux:
! 
!                               |___________S2__________|
!                               |                       |
!                       |___________S1__________|       |
!                       |                       |       |
!               |___________S0__________|       |       |
!            ...|---o---|---o---|---o---|---o---|---o---|...
!               |  i-2  |  i-1  |   i   |  i+1  |  i+2  |
!                                      -|
!                                     i+1/2
!
!   For negative flux:
!
!                       |___________S0__________|
!                       |                       |
!                       |       |___________S1__________|
!                       |       |                       |
!                       |       |       |___________S2__________|
!                    ...|---o---|---o---|---o---|---o---|---o---|...
!                       |  i-1  |   i   |  i+1  |  i+2  |  i+3  |
!                                       |+
!                                     i+1/2
!
! Li Dong <dongli@lasg.iap.ac.cn>
!
! - 2018-03-18: Initial creation.

program weno_adv_1d_case

  use netcdf

  implicit none

  real, allocatable :: x(:)         ! Cell center coordinates
  real, allocatable :: rho(:,:)     ! Tracer density being advected at cell centers
  real, allocatable :: flux_p_c(:)  ! Positive flux at cell centers
  real, allocatable :: flux_m_c(:)  ! Negative flux at cell centers
  real, allocatable :: flux_i(:)    ! Final flux at cell interfaces
  real dx                           ! Cell interval
  real :: dt = 1.0                  ! Time step size
  integer :: nx = 100               ! Cell number
  integer :: nt = 200               ! Integration time step number
  real :: u = 0.005                 ! Advection speed
  real coef                         ! dt / dx
  real, parameter :: eps = 1.0d-6   ! A small value to avoid divided-by-zero
  integer, parameter :: ns = 3      ! Stencil width

  integer i, time_step, old, new
  character(256) namelist_path
  logical is_exist

  namelist /params/ nx, nt, dt, u

  call get_command_argument(1, namelist_path)
  inquire(file=namelist_path, exist=is_exist)
  if (is_exist) then
    open(10, file=namelist_path)
    read(10, nml=params)
    close(10)
  end if

  allocate(x(nx))
  allocate(rho(1-ns:nx+ns,2))
  allocate(flux_p_c(1-ns:nx+ns))
  allocate(flux_m_c(1-ns:nx+ns))
  allocate(flux_i(1:nx+1))

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
    call weno(rho(:,old))
    do i = 1, nx
      rho(i,new) = rho(i,old) - coef * (flux_i(i+1) - flux_i(i))
    end do
    call full_boundary_condition(rho(:,new))
    ! RK 2nd stage
    call weno(rho(:,new))
    do i = 1, nx
      rho(i,new) = (3.0d0 * rho(i,old) + rho(i,new) - coef * (flux_i(i+1) - flux_i(i))) / 4.0d0
    end do
    call full_boundary_condition(rho(:,new))
    ! RK 3st stage
    call weno(rho(:,new))
    do i = 1, nx
      rho(i,new) = (rho(i,old) + 2.0d0 * (rho(i,new) - coef * (flux_i(i+1) - flux_i(i)))) / 3.0d0
    end do
    call full_boundary_condition(rho(:,new))
    ! Change time indices.
    i = old; old = new; new = i
    time_step = time_step + 1
    call output(rho(:,old))
    print *, time_step, sum(rho(1:nx,old))
  end do

  deallocate(x)
  deallocate(rho)
  deallocate(flux_p_c)
  deallocate(flux_m_c)
  deallocate(flux_i)

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

  subroutine weno(rho)

    real, intent(in) :: rho(1-ns:nx+ns)

    real, parameter :: c1 = 13.0d0 / 12.0d0
    real, parameter :: c2 = 0.25d0
    real, parameter :: gamma_p(3) = [1.0d0 / 10.0d0, 3.0d0 / 5.0d0, 3.0d0 / 10.0d0]
    real, parameter :: gamma_m(3) = [3.0d0 / 10.0d0, 3.0d0 / 5.0d0, 1.0d0 / 10.0d0]

    real alpha ! Maximum dflux/drho (e.g. advection velocity in this case)
    real beta(3), wgt(3), flux_s_i(3) ! Values on each stencil
    real flux_p_i, flux_m_i
    integer i

    alpha = abs(u) ! NOTE: In this case, max dflux/drho is just u.

    ! Calculate splitted flux at cell centers by using Lax-Friedrichs splitting.
    ! The splitting is for numerical stability by respecting upwind.
    do i = 1, nx
      flux_p_c(i) = 0.5d0 * (u + alpha) * rho(i)
      flux_m_c(i) = 0.5d0 * (u - alpha) * rho(i)
    end do
    call full_boundary_condition(flux_p_c)
    call full_boundary_condition(flux_m_c)

    do i = 1, nx
      ! Positive flux at cell interfaces
      ! - Calculate flux at interfaces for each stencil.
      flux_s_i(1) = (2 * flux_p_c(i-2) - 7 * flux_p_c(i-1) + 11 * flux_p_c(i  )) / 6.0d0
      flux_s_i(2) = (  - flux_p_c(i-1) + 5 * flux_p_c(i  ) +  2 * flux_p_c(i+1)) / 6.0d0
      flux_s_i(3) = (2 * flux_p_c(i  ) + 5 * flux_p_c(i+1)      - flux_p_c(i+2)) / 6.0d0
      ! - Calculate smooth indicators for each stencil regarding cell centers.
      beta(1) = c1 * (flux_p_c(i-2) - 2 * flux_p_c(i-1) + flux_p_c(i  ))**2 + c2 * (flux_p_c(i-2) - 4 * flux_p_c(i-1) + 3 * flux_p_c(i))**2
      beta(2) = c1 * (flux_p_c(i-1) - 2 * flux_p_c(i  ) + flux_p_c(i+1))**2 + c2 * (flux_p_c(i-1) - flux_p_c(i+1))**2
      beta(3) = c1 * (flux_p_c(i  ) - 2 * flux_p_c(i+1) + flux_p_c(i+2))**2 + c2 * (flux_p_c(i+2) - 4 * flux_p_c(i+1) + 3 * flux_p_c(i))**2
      ! - Calculate stencil linear combination weights considering smooth indicators.
      wgt(:) = gamma_p(:) / (eps + beta(:))**2
      wgt(:) = wgt(:) / sum(wgt)
      flux_p_i = sum(wgt(:) * flux_s_i(:))
      ! Negative flux at cell interfaces
      ! - Calculate flux at interfaces for each stencil.
      flux_s_i(1) = (2 * flux_m_c(i+1) + 5 * flux_m_c(i  )      - flux_m_c(i-1)) / 6.0d0
      flux_s_i(2) = (  - flux_m_c(i+2) + 5 * flux_m_c(i+1) +  2 * flux_m_c(i  )) / 6.0d0
      flux_s_i(3) = (2 * flux_m_c(i+3) - 7 * flux_m_c(i+2) + 11 * flux_m_c(i+1)) / 6.0d0
      ! - Calculate smooth indicators for each stencil regarding cell centers.
      beta(1) = c1 * (flux_m_c(i-1) - 2 * flux_m_c(i  ) + flux_m_c(i+1))**2 + c2 * (flux_m_c(i-1) - 4 * flux_m_c(i  ) + 3 * flux_m_c(i+1))**2
      beta(2) = c1 * (flux_m_c(i  ) - 2 * flux_m_c(i+1) + flux_m_c(i+2))**2 + c2 * (flux_m_c(i  ) - flux_m_c(i+2))**2
      beta(3) = c1 * (flux_m_c(i+1) - 2 * flux_m_c(i+2) + flux_m_c(i+3))**2 + c2 * (flux_m_c(i+3) - 4 * flux_m_c(i+2) + 3 * flux_m_c(i+1))**2
      ! - Calculate stencil linear combination weights considering smooth indicators.
      wgt(:) = gamma_m(:) / (eps + beta(:))**2
      wgt(:) = wgt(:) / sum(wgt)
      flux_m_i = sum(wgt(:) * flux_s_i(:))

      flux_i(i+1) = flux_p_i + flux_m_i
    end do
    call half_boundary_condition(flux_i)

  end subroutine weno

  subroutine output(rho)

    real, intent(in) :: rho(1-ns:nx+ns)

    character(30) file_name
    integer file_id, time_dim_id, time_var_id, x_dim_id, x_var_id, rho_var_id, ierr

    write(file_name, "('weno.', I3.3, '.', I4.4, '.nc')") nx, time_step

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

end program weno_adv_1d_case
