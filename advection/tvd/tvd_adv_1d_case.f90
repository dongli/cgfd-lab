! This case demonstrates the effects of TVD modifications on Lax-Wendroff scheme.
!
! Weighted-averaged flux:
!
!   f_{i+1/2} = (1 + 2 * a * c) / 2 * (u * rho_i) + (1 - 2 * a * c) / 2 * (u * rho_{i+1}).
!
! When a = 1/2, it is thhe Lax-Wendroff flux:
!
!   f_{i+1/2} = (1 + c) / 2 * (u * rho_i) + (1 - c) / 2 * (u * rho_{i+1})
!
! Li Dong <dongli@lasg.iap.ac.cn>
!
! - 2022-03-13: Initial creation.

program tvd_adv_1d_case

  use netcdf

  implicit none

  integer                           :: nx = 100   ! Cell number
  integer                           :: nt = 200   ! Integration time step number
  integer, parameter                :: ns = 1     ! Stencil width
  real, allocatable, dimension(:  ) :: x          ! Cell center coordinates
  real, allocatable, dimension(:  ) :: f          ! Flux at cell interfaces
  real, allocatable, dimension(:,:) :: rho        ! Tracer density being advected at cell centers
  real                              :: a = 0.5d0  ! Alpha parameter when it is 0.5, flux is Lax-Wendroff
  real                              :: u = 0.005d0! Advection speed
  real                              :: dx         ! Cell interval
  real                              :: dt = 1.0d0 ! Time step size
  real                              :: c          ! CFL number, u * dt / dx
  character(10) :: flux_limiter_type = 'van_leer'
  integer i
  integer :: time_step = 0, old = 1, new = 2
  character(256) namelist_path
  logical is_exist

  namelist /params/ nx, nt, dt, u, flux_limiter_type

  call get_command_argument(1, namelist_path)
  inquire(file=namelist_path, exist=is_exist)
  if (is_exist) then
    open(10, file=namelist_path)
    read(10, nml=params)
    close(10)
  end if

  allocate(x(nx))
  allocate(f(1:nx+1))
  allocate(rho(1-ns:nx+ns,2))

  ! Set mesh grid coordinates.
  dx = 1.0d0 / nx
  do i = 1, nx
    x(i) = (i - 1) * dx
  end do

  ! Set initial condition.
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
  c = u * dt / dx
  write(*, *) time_step, sum(rho(1:nx,old))
  do while (time_step < nt)
    call flux(rho(:,old), f)
    do i = 1, nx
      rho(i,new) = rho(i,old) - dt / dx * (f(i+1) - f(i))
    end do
    call full_boundary_condition(rho(:,new))
    ! Change time indices.
    i = old; old = new; new = i
    time_step = time_step + 1
    call output(rho(:,old))
    write(*, *) time_step, sum(rho(1:nx,old))
  end do

  deallocate(x, f, rho)

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

  subroutine flux(rho, f)

    real, intent(in) :: rho(1-ns:nx+ns)
    real, intent(inout) :: f(1:nx+1)

    integer i

    do i = 1, nx
      f(i+1) = u * rho(i) + 0.5d0 * u * (1 - c) * (rho(i+1) - rho(i)) * flux_limiter(rho(i-1:i+1))
    end do
    call half_boundary_condition(f)

  end subroutine flux

  real function flux_limiter(rho) result(res)

    real, intent(in) :: rho(-1:1)

    real, parameter :: eps = 1.0e-6
    real r

    r = (rho(0) - rho(-1) + eps) / (rho(1) - rho(0) + eps)

    select case (flux_limiter_type)
    case ('none')
      res = 1
    case ('minmod')
      res = max(0.0, min(1.0, r))
    case ('superbee')
      res = max(0.0, min(1.0, 2 * r), min(2.0, r))
    case ('van_leer')
      res = (r + abs(r)) / (1 + abs(r))
    case ('mc')
      res = max(0.0, min(2 * r, (1 + r) / 2, 2.0))
    case default
      write(*, *) '[Error]: Invalide flux_limiter_type ' // trim(flux_limiter_type) // '!'
    end select

  end function flux_limiter

  subroutine output(rho)

    real, intent(in) :: rho(1-ns:nx+ns)

    character(30) file_name
    integer ncid, time_dimid, time_varid, x_dimid, x_varid, rho_varid, ierr

    write(file_name, "('tvd.', A, '.', I3.3, '.', I4.4, '.nc')") trim(flux_limiter_type), nx, time_step

    ierr = NF90_CREATE(file_name, NF90_CLOBBER, ncid)
    ierr = NF90_PUT_ATT(ncid, NF90_GLOBAL, 'scheme', 'Lax-Wendroff (' // trim(flux_limiter_type) // ')')

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

end program tvd_adv_1d_case
