module linear_swm_common_mod

  use iso_fortran_env, only: rp => real64

  implicit none

  integer, parameter :: nx = 10, ny = 10, hw = 1
  real(rp), allocatable :: full_x(:), half_x(:)
  real(rp), allocatable :: full_y(:), half_y(:)
  real(rp), allocatable :: u(:,:,:), v(:,:,:), h(:,:,:)
  real(rp), allocatable :: du(:,:,:), dv(:,:,:), dh(:,:,:)

  real(rp), parameter :: pi = atan(1.0_rp) * 4.0_rp
  real(rp), parameter :: omg = 2.0_rp * pi / 86400.0_rp
  real(rp), parameter :: f = 2.0_rp * omg * sin(pi * 0.25_rp)
  real(rp), parameter :: g = 9.8_rp
  real(rp), parameter :: h0 = 0.1_rp

  real(rp) dx, dy, dt
  integer i, j, time_step, old, new

contains

  subroutine init_swm()

    integer i, j

    allocate(full_x(1-hw:nx+hw), half_x(1-hw:nx+hw))
    allocate(full_y(1-hw:ny+hw), half_y(1-hw:ny+hw))
    allocate(u(1-hw:nx+hw,1-hw:ny+hw,0:2))
    allocate(v(1-hw:nx+hw,1-hw:ny+hw,0:2))
    allocate(h(1-hw:nx+hw,1-hw:ny+hw,0:2))
    allocate(du(1-hw:nx+hw,1-hw:ny+hw,0:2))
    allocate(dv(1-hw:nx+hw,1-hw:ny+hw,0:2))
    allocate(dh(1-hw:nx+hw,1-hw:ny+hw,0:2))

    time_step = 0
    old = 1
    new = 2

    ! Set mesh grids.
    dx = 1.0_rp / nx
    do i = lbound(full_x, 1), ubound(full_x, 1)
      full_x(i) = (i - 1) * dx
      half_x(i) = (i - 0.5_rp) * dx
    end do
    dy = 1.0_rp / ny
    do j = lbound(full_y, 1), ubound(full_y, 1)
      full_y(j) = (j - 1) * dy
      half_y(j) = (j - 0.5_rp) * dy
    end do

  end subroutine init_swm

  subroutine set_ic_h_spike(k)

    integer, intent(in) :: k

    u(:,:,k) = 0.0_rp
    v(:,:,k) = 0.0_rp
    h(:,:,k) = 0.0_rp
    h(nx/2,ny/2,k) = 0.1_rp * h0
    call apply_bc(k)

  end subroutine set_ic_h_spike

  subroutine apply_bc(k)

    integer, intent(in) :: k

    integer i, j

    ! u
    do i = lbound(u, 1), 0
      u(i,:,k) = u(nx+i,:,k)
    end do
    do i = nx + 1, ubound(u, 1)
      u(i,:,k) = u(i-nx,:,k)
    end do
    do j = lbound(u, 2), 0
      u(:,j,k) = u(:,ny+j,k)
    end do
    do j = ny + 1, ubound(u, 2)
      u(:,j,k) = u(:,j-ny,k)
    end do
    ! v
    do i = lbound(v, 1), 0
      v(i,:,k) = v(nx+i,:,k)
    end do
    do i = nx + 1, ubound(v, 1)
      v(i,:,k) = v(i-nx,:,k)
    end do
    do j = lbound(v, 2), 0
      v(:,j,k) = v(:,ny+j,k)
    end do
    do j = ny + 1, ubound(v, 2)
      v(:,j,k) = v(:,j-ny,k)
    end do
    ! h
    do i = lbound(h, 1), 0
      h(i,:,k) = h(nx+i,:,k)
    end do
    do i = nx + 1, ubound(h, 1)
      h(i,:,k) = h(i-nx,:,k)
    end do
    do j = lbound(h, 2), 0
      h(:,j,k) = h(:,ny+j,k)
    end do
    do j = ny + 1, ubound(h, 2)
      h(:,j,k) = h(:,j-ny,k)
    end do

  end subroutine apply_bc

  subroutine time_advance()

    integer tmp

    tmp = old
    old = new
    new = tmp

    time_step = time_step + 1

  end subroutine time_advance

end module linear_swm_common_mod
