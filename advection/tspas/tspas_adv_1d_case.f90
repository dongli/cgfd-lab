! This is a 1D advection example using square initial condition and periodic
! boundary condition for TSPAS finite difference scheme.
!
! Li Dong <dongli@lasg.iap.ac.cn>
!
! - 2018-03-21: Initial creation.

program tspas_adv_1d_case

  use netcdf

  implicit none

  real, allocatable :: x(:)         ! Cell center coordinates
  real, allocatable :: xi(:)        ! Cell interface coordinates
  real, allocatable :: rho(:,:)     ! Tracer density being advected at cell centers
  real, allocatable :: flux(:)      ! Flux at cell interfaces
  real, allocatable :: gamma(:)     !
  real, allocatable :: A(:)         ! 
  real, allocatable :: c_star(:)    ! Switch of upwind scheme
  real, allocatable :: u_star(:)    ! Modified velocity
  real dx                           ! Cell interval
  real :: dt = 1.0                  ! Time step size
  integer :: nx = 100               ! Cell number
  integer :: nt = 200               ! Integration time step number
  real :: u = 0.005                 ! Advection speed
  real coef                         ! dt / dx
  real, parameter :: eps = 1.0d-80  ! A small value to avoid divided-by-zero
  integer, parameter :: ns = 1      ! Stencil width
  integer, parameter :: star = 0
  integer i
  integer :: time_step = 0, old = 1, new = 2
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
  allocate(xi(nx+1))
  allocate(rho(1-ns:nx+ns,0:2))
  allocate(flux(1:nx+1))
  allocate(gamma(1:nx+1))
  allocate(A(1-ns:nx+ns))
  allocate(c_star(1:nx+1))
  allocate(u_star(1:nx+1))

  ! Set mesh grid coordinates.
  dx = 1.0d0 / nx
  do i = 1, nx
    x(i) = (i - 1) * dx
    xi(i) = x(i) - 0.5d0 * dx
  end do
  xi(nx+1) = x(i) + 0.5d0 * dx

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
    call tspas(rho(:,old), rho(:,star))
    do i = 1, nx
      rho(i,new) = rho(i,old) - coef * (flux(i+1) - flux(i))
    end do
    call full_boundary_condition(rho(:,new))
    ! Change time indices.
    i = old; old = new; new = i
    time_step = time_step + 1
    call output(rho(:,old))
    print *, time_step, sum(rho(1:nx,old))
  end do

  deallocate(x)
  deallocate(xi)
  deallocate(rho)
  deallocate(flux)
  deallocate(gamma)
  deallocate(A)
  deallocate(c_star)
  deallocate(u_star)

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

  subroutine tspas(rho, rho_star)

    real, intent(in) :: rho(1-ns:nx+ns)
    real, intent(inout) :: rho_star(1-ns:nx+ns)

    real alpha, beta, tmp1, tmp2, tmp3
    integer i

    ! Run Lax-Wendroff pass.
    do i = 1, nx
      flux(i+1) = 0.5d0 * (u * (rho(i+1) + rho(i)) - coef * u**2 * (rho(i+1) - rho(i)))
    end do
    call half_boundary_condition(flux)
    do i = 1, nx + 1
      alpha = abs(u) * coef
      gamma(i) = alpha * (1.0d0 - alpha)
    end do
    do i = 1, nx
      beta = max(2.0d0 / (2.0d0 - gamma(i-1)), 2.0d0 / (2.0d0 - gamma(i)))
      rho_star(i) = rho(i) - beta * coef * (flux(i+1) - flux(i))
    end do
    call full_boundary_condition(rho_star)
    ! Calculate A.
    do i = 1, nx
      A(i) = (rho_star(i) - max(rho(i-1), rho(i), rho(i+1))) * &
             (rho_star(i) - min(rho(i-1), rho(i), rho(i+1)))
    end do
    call full_boundary_condition(A)
    ! Calculate u_star
    do i = 1, nx
      tmp1 = abs(A(i  )) + A(i  )
      tmp2 = abs(A(i+1)) + A(i+1)
      tmp3 = 0.5d0  * (tmp1 / (abs(A(i)) + eps) + tmp2 / (abs(A(i+1)) + eps)) - &
             0.25d0 * (tmp1 * tmp2 / (abs(A(i) * A(i+1)) + eps))
      c_star(i+1) = tmp3
      u_star(i+1) = (c_star(i+1) + (1 - c_star(i+1)) * coef * abs(u)) * u
    end do
    ! Run upwind pass.
    do i = 1, nx
      flux(i+1) = 0.5d0 * (u * (rho(i+1) + rho(i)) - abs(u_star(i+1)) * (rho(i+1) - rho(i)))
    end do
    call half_boundary_condition(flux)

  end subroutine tspas

  subroutine output(rho)

    real, intent(in) :: rho(1-ns:nx+ns)

    character(30) file_name
    integer ncid, time_dimid, time_varid, x_dimid, x_varid
    integer xi_dimid, xi_varid, rho_varid, c_star_varid, ierr

    write(file_name, "('tspas.', I3.3, '.', I4.4, '.nc')") nx, time_step

    ierr = NF90_CREATE(file_name, NF90_CLOBBER, ncid)
    ierr = NF90_PUT_ATT(ncid, NF90_GLOBAL, 'scheme', 'TSPAS')

    ierr = NF90_DEF_DIM(ncid, 'time', NF90_UNLIMITED, time_dimid)
    ierr = NF90_DEF_VAR(ncid, 'time', NF90_INT, [time_dimid], time_varid)
    ierr = NF90_DEF_DIM(ncid, 'x', nx, x_dimid)
    ierr = NF90_DEF_VAR(ncid, 'x', NF90_FLOAT, [x_dimid], x_varid)
    ierr = NF90_DEF_DIM(ncid, 'xi', nx + 1, xi_dimid)
    ierr = NF90_DEF_VAR(ncid, 'xi', NF90_FLOAT, [xi_dimid], xi_varid)
    ierr = NF90_DEF_VAR(ncid, 'rho', NF90_FLOAT, [x_dimid, time_dimid], rho_varid)
    ierr = NF90_DEF_VAR(ncid, 'c_star', NF90_FLOAT, [xi_dimid, time_dimid], c_star_varid)

    ierr = NF90_ENDDEF(ncid)

    ierr = NF90_PUT_VAR(ncid, time_varid, time_step)
    ierr = NF90_PUT_VAR(ncid, x_varid, x)
    ierr = NF90_PUT_VAR(ncid, xi_varid, xi)
    ierr = NF90_PUT_VAR(ncid, rho_varid, rho(1:nx))
    ierr = NF90_PUT_VAR(ncid, c_star_varid, c_star(1:nx+1))

    ierr = NF90_CLOSE(ncid)

  end subroutine output

end program tspas_adv_1d_case
