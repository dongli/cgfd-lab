! This is a 1D advection example using square initial condition and periodic
! boundary condition for Crank-Nicolson finite difference scheme.
!
! Li Dong <dongli@lasg.iap.ac.cn>
!
! - 2018-03-20: Initial creation.

program crank_nicolson_adv_1d_case

  use netcdf
  use iso_c_binding

  implicit none

  interface
    type(c_ptr) function gsl_vector_alloc(n) bind(c) result(vec)
      use iso_c_binding
      integer(c_size_t), value :: n
    end function gsl_vector_alloc

    subroutine gsl_vector_set(vec, idx, val) bind(c)
      use iso_c_binding
      type(c_ptr), value :: vec
      integer(c_size_t), value :: idx
      real(c_double), value :: val
    end subroutine gsl_vector_set

    real(c_double) function gsl_vector_get(vec, idx) bind(c)
      use iso_c_binding
      type(c_ptr), value :: vec
      integer(c_size_t), value :: idx
    end function gsl_vector_get

    subroutine gsl_vector_free(vec) bind(c)
      use iso_c_binding
      type(c_ptr), value :: vec
    end subroutine gsl_vector_free

    subroutine gsl_linalg_solve_cyc_tridiag(diag, e, f, b, x) bind(c)
      use iso_c_binding
      type(c_ptr), value :: diag
      type(c_ptr), value :: e
      type(c_ptr), value :: f
      type(c_ptr), value :: b
      type(c_ptr), value :: x
    end subroutine gsl_linalg_solve_cyc_tridiag
  end interface

  real, allocatable :: x(:)         ! Cell center coordinates
  real, allocatable :: rho(:,:)     ! Tracer density being advected at cell centers

  type(c_ptr) a1
  type(c_ptr) a2
  type(c_ptr) a3
  type(c_ptr) b
  type(c_ptr) y

  real dx                           ! Cell interval
  real :: dt = 1.0                  ! Time step size
  integer :: nx = 100               ! Cell number
  integer :: nt = 200               ! Integration time step number
  real :: u = 0.005                 ! Advection speed
  real coef                         ! dt / dx * 0.25d0
  integer, parameter :: ns = 1      ! Stencil width
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
  allocate(rho(1-ns:nx+ns,2))
  a1 = gsl_vector_alloc(int(nx, c_size_t))
  a2 = gsl_vector_alloc(int(nx, c_size_t))
  a3 = gsl_vector_alloc(int(nx, c_size_t))
  b = gsl_vector_alloc(int(nx, c_size_t))
  y = gsl_vector_alloc(int(nx, c_size_t))

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
  coef = dt / dx * 0.25d0
  print *, time_step, sum(rho(1:nx,old))
  do while (time_step < nt)
    call crank_nicolson(rho(:,old))
    do i = 1, nx
      rho(i,new) = gsl_vector_get(y, int(i - 1, c_size_t))
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
  call gsl_vector_free(a1)
  call gsl_vector_free(a2)
  call gsl_vector_free(a3)
  call gsl_vector_free(b)
  call gsl_vector_free(y)

contains

  subroutine full_boundary_condition(x)

    real, intent(inout) :: x(1-ns:nx+ns)

    integer i

    do i = 1, ns
      x(1-i) = x(nx-i+1)
      x(nx+i) = x(1+i-1)
    end do

  end subroutine full_boundary_condition

  subroutine crank_nicolson(rho)

    real, intent(in) :: rho(1-ns:nx+ns)

    integer(c_size_t) i

    do i = 1, nx
      ! Although u is constant in this case, we still write as it is variant.
      if (i == 1) then
        call gsl_vector_set(a1, i - 1, - coef * u)
        call gsl_vector_set(a3, i - 1,   coef * u)
      else if (i == nx) then
        call gsl_vector_set(a1, i - 1, - coef * u)
        call gsl_vector_set(a3, i - 1,   coef * u)
      else
        call gsl_vector_set(a1, i - 1, - coef * u)
        call gsl_vector_set(a3, i - 1,   coef * u)
      end if
      call gsl_vector_set(a2, i - 1, 1.0)
      call gsl_vector_set(b, i - 1, coef * u * rho(i-1) + rho(i) - coef * u * rho(i+1))
    end do
    call gsl_linalg_solve_cyc_tridiag(a2, a3, a1, b, y)

  end subroutine crank_nicolson

  subroutine output(rho)

    real, intent(in) :: rho(1-ns:nx+ns)

    character(30) file_name
    integer ncid, time_dimid, time_varid, x_dimid, x_varid, rho_varid, ierr

    write(file_name, "('crank_nicolson.', I3.3, '.', I4.4, '.nc')") nx, time_step

    ierr = NF90_CREATE(file_name, NF90_CLOBBER, ncid)
    ierr = NF90_PUT_ATT(ncid, NF90_GLOBAL, 'scheme', 'Crank-Nicolson')

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

end program crank_nicolson_adv_1d_case
