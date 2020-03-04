! This is a 1D advection example using square initial condition and periodic
! boundary condition for MPDATA finite difference scheme.
!
! Li Dong <dongli@lasg.iap.ac.cn>
!
! - 2018-03-21: Initial creation.

program mpdata_adv_1d_case

  use netcdf

  implicit none

  real, allocatable :: x(:)         ! Cell center coordinates
  real, allocatable :: rho(:,:)     ! Tracer density being advected at cell centers
  real, allocatable :: flux(:)      ! Flux at cell interfaces
  real, allocatable :: uc(:)        ! Antidiffusion velocity
  real dx                           ! Cell interval
  real :: dt = 1.0                  ! Time step size
  integer :: nx = 100               ! Cell number
  integer :: nt = 200               ! Integration time step number
  integer :: iord = 3               ! Scheme order
  real :: u = 0.005                 ! Advection speed
  real coef                         ! dt / dx
  real, parameter :: eps = 1.0d-15  ! A small value to avoid divided-by-zero
  integer, parameter :: ns = 1      ! Stencil width
  integer i
  integer :: time_step = 0, old = 1, new = 2
  character(256) namelist_path
  logical is_exist

  namelist /params/ nx, nt, dt, iord, u

  call get_command_argument(1, namelist_path)
  inquire(file=namelist_path, exist=is_exist)
  if (is_exist) then
    open(10, file=namelist_path)
    read(10, nml=params)
    close(10)
  end if

  allocate(x(nx))
  allocate(rho(1-ns:nx+ns,old:new))
  allocate(flux(1:nx+1))
  allocate(uc(1:nx+1))

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
  coef = dt / dx
  print *, time_step, sum(rho(1:nx,old))
  do while (time_step < nt)
    call mpdata(rho(:,old), rho(:,new))
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

  subroutine mpdata(rho, rho_star)

    real, intent(in) :: rho(1-ns:nx+ns)
    real, intent(inout) :: rho_star(1-ns:nx+ns)

    integer i, j

    uc(:) = u
    rho_star(:) = rho(:)
    do j = 1, iord
      if (j > 1) then
        ! Calculate antidiffusion velocity.
        do i = 1, nx
          uc(i+1) = (abs(uc(i+1)) - uc(i+1)**2 * coef) * &
                    (rho_star(i+1) - rho_star(i)) / (rho_star(i+1) + rho_star(i) + eps)
        end do
      end if
      call upwind(rho_star, uc)
      ! Update tracer density.
      do i = 1, nx
        rho_star(i) = rho_star(i) - coef * (flux(i+1) - flux(i))
      end do
      call full_boundary_condition(rho_star)
    end do

  end subroutine mpdata

  subroutine upwind(rho, u)

    real, intent(in) :: rho(1-ns:nx+ns)
    real, intent(in) :: u(1:nx+1)

    integer i

    do i = 1, nx
      flux(i+1) = 0.5d0 * (u(i+1) * (rho(i+1) + rho(i)) - abs(u(i+1)) * (rho(i+1) - rho(i)))
    end do
    call half_boundary_condition(flux)

  end subroutine upwind

  subroutine output(rho)

    real, intent(in) :: rho(1-ns:nx+ns)

    character(30) file_name
    integer ncid, time_dimid, time_varid, x_dimid, x_varid, rho_varid, ierr

    write(file_name, "('mpdata.', I3.3, '.', I4.4, '.nc')") nx, time_step

    ierr = NF90_CREATE(file_name, NF90_CLOBBER, ncid)
    ierr = NF90_PUT_ATT(ncid, NF90_GLOBAL, 'scheme', 'MPDATA')

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

end program mpdata_adv_1d_case
