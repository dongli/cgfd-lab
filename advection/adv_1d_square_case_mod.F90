module adv_1d_square_case_mod

  use adv_1d_common_mod

  implicit none

  real, allocatable :: x(:)        ! Cell center coordinates
  real dx                          ! Cell interval
  real :: dt = 1.0                 ! Time step size
  integer :: nx = 100              ! Cell number
  integer :: nt = 200              ! Integration time step number
  real :: u = 0.005                ! Advection speed

contains

  subroutine adv_1d_square_case_init(ns, rho)

    integer, intent(in) :: ns
    real, intent(out) :: rho(1-ns+1:nx+ns)

    integer i

    allocate(x(nx))

    ! Set mesh grid coordinates.
    dx = 1.0d0 / nx
    do i = 1, nx
      x(i) = (i - 1) * dx
    end do

    ! Set initial condition.
    do i = 1, nx
      if (x(i) >= 0.05 .and. x(i) <= 0.3) then
        rho(i) = 1.0d0
      else
        rho(i) = 0.0d0
      end if
    end do
    call apply_bc(ns, nx, rho)

  end subroutine adv_1d_square_case_init

  subroutine adv_1d_square_case_final()

    deallocate(x)

  end subroutine adv_1d_square_case_final

end module adv_1d_square_case_mod