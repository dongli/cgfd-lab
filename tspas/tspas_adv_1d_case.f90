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
  real, allocatable :: rho(:,:)     ! Tracer density being advected at cell centers
  real, allocatable :: flux(:)      ! Flux at cell interfaces
  real, allocatable :: gamma(:)     !
  real, allocatable :: A(:)         ! 
  real, allocatable :: c_star(:)    ! Switch of upwind scheme
  real, allocatable :: u_star(:)    ! Modified velocity
  real dx                           ! Cell interval
  real dt                           ! Time step size
  integer nx                        ! Cell number
  integer nt                        ! Integration time step number

  real :: u = 0.005                 ! Advection speed
  real coef                         ! dt / dx
  real, parameter :: eps = 1.0d-80  ! A small value to avoid divided-by-zero
  integer, parameter :: ns = 1      ! Stencil width
  integer, parameter :: star = 0

  integer i, time_step, old, new
  character(256) namelist_path
  logical is_exist

  namelist /params/ nx, nt, dx, dt, u

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
    integer file_id, time_dim_id, time_var_id, x_dim_id, x_var_id, rho_var_id, ierr

    write(file_name, "('tspas.', I3.3, '.nc')") time_step

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

end program tspas_adv_1d_case
