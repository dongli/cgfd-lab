! This is a 2D advection example using slotted cylinder initial condition and
! double periodic boundary condition for semi-Lagrangian scheme.
!
! Li Dong <dongli@lasg.iap.ac.cn>
!
! - 2018-03-25: Initial creation.

program semi_lagrangian_adv_2d_case

  use netcdf

  implicit none

  real, allocatable :: x(:)               ! Cell center coordinates
  real, allocatable :: y(:)               ! Cell center coordinates
  real, allocatable :: u(:,:)             ! Velocity component along x axis
  real, allocatable :: v(:,:)             ! Velocity component along y axis
  real, allocatable :: rho(:,:,:)         ! Tracer density being advected at cell centers
  real dx                                 ! Cell interval along x axis
  real dy                                 ! Cell interval along y axis
  real :: dt = 0.8                        ! Time step size
  integer :: nx = 100                     ! Cell number along x axis
  integer :: ny = 100                     ! Cell number along y axis
  integer :: nt = 263                     ! Integration time step number
  real, parameter :: omega = 0.03         ! Rotation angular speed
  real, parameter :: x0 = 0.25            ! Initial coordinate
  real, parameter :: y0 = 0.5             ! Initial coordinate
  integer, parameter :: ns = 2            ! Stencil width
  integer i, j
  integer :: time_step = 0, old = 1, new = 2
  integer :: velocity_interp_type = 1     ! 1: bilinear; 2: biquadratic; 3: bicubic
  integer :: density_interp_type = 3      ! 1: bilinear; 2: biquadratic; 3: bicubic
  character(256) namelist_path
  logical is_exist

  namelist /params/ nx, dx, ny, dy, nt, dt, velocity_interp_type, density_interp_type

  call get_command_argument(1, namelist_path)
  inquire(file=namelist_path, exist=is_exist)
  if (is_exist) then
    open(10, file=namelist_path)
    read(10, nml=params)
    close(10)
  end if

  allocate(x(1-ns:nx+ns))
  allocate(y(1-ns:ny+ns))
  allocate(u(1-ns:nx+ns,1-ns:ny+ns))
  allocate(v(1-ns:nx+ns,1-ns:ny+ns))
  allocate(rho(1-ns:nx+ns,1-ns:ny+ns,old:new))

  ! Set mesh grid coordinates.
  dx = 1.0d0 / nx
  do i = 1, nx
    x(i) = (i - 1) * dx
  end do
  do i = 1, ns
    x(1-i) = x(1) - i * dx
    x(nx+i) = x(nx) + i * dx
  end do
  dy = 1.0d0 / ny
  do j = 1, ny
    y(j) = (j - 1) * dy
  end do
  do j = 1, ns
    y(1-j) = y(1) - j * dy
    y(ny+j) = y(ny) + j * dy
  end do

  ! Set initial condition.
  do j = 1, ny
    do i = 1, nx
      u(i,j) = - omega * (y(j) - 0.5)
      v(i,j) =   omega * (x(i) - 0.5)
      ! u(i,j) = 0.005
      ! v(i,j) = 0.005
      if (abs(x(i) - x0) >= 0.02 .and. (x(i) - x0)**2 + (y(j) - y0)**2 < 0.01) then
        rho(i,j,old) = 1.0
      else if (abs(x(i) - x0) < 0.02 .and. (x(i) - x0)**2 + (y(j) - y0)**2 < 0.01 .and. y(j) >= 0.55) then
        rho(i,j,old) = 1.0
      else
        rho(i,j,old) = 0.0
      end if
    end do
  end do
  call full_boundary_condition(u)
  call full_boundary_condition(v)
  call full_boundary_condition(rho(:,:,old))
  call output(rho(:,:,old))

  ! Run integration.
  print *, time_step, sum(rho(1:nx,1:ny,old))
  do while (time_step < nt)
    call semi_lagrangian()
    ! Change time indices.
    i = old; old = new; new = i
    time_step = time_step + 1
    call output(rho(:,:,old))
    print *, time_step, sum(rho(1:nx,1:ny,old))
  end do

  deallocate(x)
  deallocate(y)
  deallocate(u)
  deallocate(v)
  deallocate(rho)

contains

  subroutine full_boundary_condition(array)

    real, intent(inout) :: array(1-ns:nx+ns,1-ns:ny+ns)

    integer i, j

    do j = 1, ny
      do i = 1, ns
        array(1-i,j) = array(nx-i+1,j)
        array(nx+i,j) = array(i,j)
      end do
    end do
    do j = 1, ns
      do i = 1, nx
        array(i,1-j) = array(i,ny-j+1)
        array(i,ny+j) = array(i,j)
      end do
    end do

  end subroutine full_boundary_condition

  subroutine semi_lagrangian()

    real um, vm, xd, yd, dx, dy, xm, ym
    integer i, j, iter

    do j = 1, ny
      do i = 1, nx
        xd = x(i)
        yd = y(j)
        do iter = 1, 2
          call diff_coord(x(i), y(j), xd, yd, dx, dy)
          call add_coord(x(i), y(j), -0.5 * dx, -0.5 * dy, xm, ym)
          call interp_velocity(u, v, xm, ym, um, vm)
          call add_coord(x(i), y(j), -dt * um, -dt * vm, xd, yd)
        end do
        call interp_density(rho(:,:,old), xd, yd, rho(i,j,new))
      end do
    end do
    call full_boundary_condition(rho(:,:,new))

  end subroutine semi_lagrangian

  subroutine diff_coord(x1, y1, x2, y2, dx, dy)

    real, intent(in) :: x1, y1, x2, y2
    real, intent(out) :: dx, dy

    real dx1, dy1, dx2, dy2

    dx = x1 - x2
    dy = y1 - y2
    dx1 = abs(dx)
    dy1 = abs(dy)
    dx2 = 1.0 - dx1
    dy2 = 1.0 - dy1
    if (dx1 > dx2) dx = -sign(dx2, dx)
    if (dy1 > dy2) dy = -sign(dy2, dy)

  end subroutine diff_coord

  subroutine add_coord(x1, y1, dx, dy, x2, y2)

    real, intent(in) :: x1, y1, dx, dy
    real, intent(out) :: x2, y2

    x2 = x1 + dx
    y2 = y1 + dy

    if (x2 < 0.0) x2 = 1.0 + x2
    if (x2 > 1.0) x2 = x2 - 1.0
    if (y2 < 0.0) y2 = 1.0 + y2
    if (y2 > 1.0) y2 = y2 - 1.0

  end subroutine add_coord

  subroutine interp_velocity(u, v, x0, y0, u0, v0)

    real, intent(in) :: u(1-ns:nx+ns,1-ns:ny+ns)
    real, intent(in) :: v(1-ns:nx+ns,1-ns:ny+ns)
    real, intent(in) :: x0, y0
    real, intent(out) :: u0, v0

    call interp(velocity_interp_type, u, x0, y0, u0)
    call interp(velocity_interp_type, v, x0, y0, v0)

  end subroutine interp_velocity

  subroutine interp_density(rho, x0, y0, rho0)

    real, intent(in) :: rho(1-ns:nx+ns,1-ns:ny+ns)
    real, intent(in) :: x0, y0
    real, intent(out) :: rho0

    call interp(density_interp_type, rho, x0, y0, rho0)

  end subroutine interp_density

  subroutine interp(method, f, x0, y0, f0)

    integer, intent(in) :: method
    real, intent(in) :: f(1-ns:nx+ns,1-ns:ny+ns)
    real, intent(in) :: x0, y0
    real, intent(out) :: f0

    integer n, i(4), j(4), p, q
    real wx(4), wy(4), x1, y1, x2, y2

    select case (method)
    case (1)
      n = 2
    case (2)
      n = 3
    case (3)
      n = 4
    end select
    ! Set indices of involved grids.
    i(1) = int(x0 / dx) + 1 - n / 2 + 1
    j(1) = int(y0 / dy) + 1 - n / 2 + 1
    do p = 1, n - 1
      i(p+1) = i(p) + 1
      j(p+1) = j(p) + 1
    end do
    ! Calculate weights.
    do p = 1, n
      x1 = x(i(p))
      y1 = y(j(p))
      wx(p) = 1.0
      wy(p) = 1.0
      do q = 1, n
        if (p == q) cycle
        x2 = x(i(q))
        y2 = y(j(q))
        wx(p) = wx(p) * (x0 - x2) / (x1 - x2)
        wy(p) = wy(p) * (y0 - y2) / (y1 - y2)
      end do
    end do
    ! Evaluate formula.
    f0 = 0.0
    do p = 1, n
      do q = 1, n
        f0 = f0 + wx(p) * wy(q) * f(i(p), j(q))
      end do
    end do

  end subroutine interp

  subroutine output(rho)

    real, intent(in) :: rho(1-ns:nx+ns,1-ns:ny+ns)

    character(30) file_name
    integer file_id, ierr
    integer time_dim_id, time_var_id
    integer x_dim_id, x_var_id
    integer y_dim_id, y_var_id
    integer rho_var_id

    write(file_name, "('semi_lagrangian.', I3.3, '.nc')") time_step

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

    ierr = nf90_def_dim(file_id, 'y', nx, y_dim_id)
    if (ierr /= nf90_noerr) then
      write(*, *) '[Error]: Failed to define y dimension!'
      stop 1
    end if

    ierr = nf90_def_var(file_id, 'y', nf90_float, [y_dim_id], y_var_id)
    if (ierr /= nf90_noerr) then
      write(*, *) '[Error]: Failed to define y variable!'
      stop 1
    end if

    ierr = nf90_def_var(file_id, 'rho', nf90_float, [x_dim_id, y_dim_id, time_dim_id], rho_var_id)
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

    ierr = nf90_put_var(file_id, x_var_id, x(1:nx))
    if (ierr /= nf90_noerr) then
      write(*, *) '[Error]: Failed to write x variable!'
      stop 1
    end if

    ierr = nf90_put_var(file_id, y_var_id, y(1:ny))
    if (ierr /= nf90_noerr) then
      write(*, *) '[Error]: Failed to write y variable!'
      stop 1
    end if

    ierr = nf90_put_var(file_id, rho_var_id, rho(1:nx,1:ny))
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

end program semi_lagrangian_adv_2d_case