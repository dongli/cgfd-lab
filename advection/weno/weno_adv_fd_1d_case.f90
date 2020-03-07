! This is a 1D advection example using square initial condition and periodic
! boundary condition for WENO finite difference scheme.
!
! Stencils:
!
!   For positive flux:
! 
!                               |-----------S3----------|
!                       |-----------S2----------|       |
!               |-----------S1----------|       |       |
!            ...|---o---|---o---|---o---|---o---|---o---|...
!               |  i-2  |  i-1  |   i   |  i+1  |  i+2  |
!                                      -|
!                                     i+1/2
!
!   For negative flux:
!
!                       |-----------S3----------|
!                       |       |-----------S2----------|
!                       |       |       |-----------S1----------|
!                    ...|---o---|---o---|---o---|---o---|---o---|...
!                       |  i-1  |   i   |  i+1  |  i+2  |  i+3  |
!                                       |+
!                                     i+1/2
!
! First, define a function h(x) from flux f(x):
!
!   _      1   / x_{i+1/2}
!   h_i = ---  |           h(𝜉) d𝜉 = f(x_i)
!         d x  |
!              / x_{i-1/2}
!
! which is known. Take derivative respect to x on both sides
!
!    1                                  df |
!   --- (h(x_{i+1/2}) - h(x_{i-1/2})) = -- |
!   Δ x                                 dx | x=x_i
!
! then numerical flux at i+1/2 cell interface is
!
!   ^
!   f_{i+1/2} = h(x_{i+1/2})
!
! The Lax-Friedrichs flux splitting technique is used to accommodate difference
! signs of flow velocity
!
!        +    -
!   f = f  + f
!
! where
!
!    +    1                   -    1
!   f  = --- (f + umax rho), f  = --- (f - umax rho)
!         2                        2
!
! References:
!
! - Chi-Wang Shu, 2009: High Order Weighted Essentially Non-Oscillatory Schemes for Convection Dominated Problems. SIAM Rev., 51(1).
!
! Li Dong <dongli@lasg.iap.ac.cn>
!
! - 2018-03-18: Initial creation.

program weno_adv_fd_1d_case

  use adv_1d_square_case_mod

  implicit none

  real, allocatable, dimension(:,:) :: rho  ! Tracer density being advected at cell centers
  real, allocatable, dimension(:,:) :: dfdx ! Tendency at cell centers
  real, allocatable, dimension(:  ) :: fpc  ! Positive flux at cell centers
  real, allocatable, dimension(:  ) :: fmc  ! Negative flux at cell centers
  real, allocatable, dimension(:  ) :: f    ! Final flux at cell interfaces
  
  real, parameter :: eps = 1.0d-16 ! A small value to avoid divided-by-zero
  integer, parameter :: ns = 3     ! Stencil width


  allocate(rho (1-ns:nx+ns,2))
  allocate(dfdx(1   :nx   ,4))
  allocate(fpc (1-ns:nx+ns  ))
  allocate(fmc (1-ns:nx+ns  ))
  allocate(f   (0   :nx+1   ))

  call adv_1d_square_case_init(ns, rho(:,old))
  call output('weno', time_step, ns, nx, x, rho(:,old))

  ! Run integration.
  write(*, *) time_step, sum(rho(1:nx,old))
  do while (time_step < nt)
    call rk4()
    ! call rk3_tvd()
    call advance_time()
    call output('weno', time_step, ns, nx, x, rho(:,old))
    write(*, *) time_step, sum(rho(1:nx,old))
  end do

  deallocate(rho)
  deallocate(dfdx)
  deallocate(fpc)
  deallocate(fmc)
  deallocate(f)

  call adv_1d_square_case_final()

contains

  subroutine weno(rho, dfdx)

    real, intent(in) :: rho(1-ns:nx+ns)
    real, intent(out) :: dfdx(1:nx)

    real, parameter :: p1(3) = [ 1.0d0 / 3.0d0, -7.0d0 / 6.0d0, 11.0d0 / 6.0d0]
    real, parameter :: p2(3) = [-1.0d0 / 6.0d0,  5.0d0 / 6.0d0,  1.0d0 / 3.0d0]
    real, parameter :: p3(3) = [ 1.0d0 / 3.0d0,  5.0d0 / 6.0d0, -1.0d0 / 6.0d0]
    real, parameter :: g (3) = [ 1.0d0 / 10.0d0, 3.0d0 / 5.0d0, 3.0d0 / 10.0d0]
    real, parameter :: s1 = 13.0d0 / 12.0d0
    real, parameter :: s2 = 0.25d0
    real, parameter :: s11(3) = [ 1.0d0, -2.0d0,  1.0d0]
    real, parameter :: s12(3) = [ 1.0d0, -4.0d0,  3.0d0]
    real, parameter :: s21(3) = [ 1.0d0, -2.0d0,  1.0d0]
    real, parameter :: s22(3) = [ 1.0d0,  0.0d0, -1.0d0]
    real, parameter :: s31(3) = [ 1.0d0, -2.0d0,  1.0d0]
    real, parameter :: s32(3) = [ 3.0d0, -4.0d0,  1.0d0]

    real umax  ! Maximum df/drho (e.g. advection velocity in this case)
    real b(3)  ! Smooth indicators
    real w(3)  ! Combination weights
    real fs(3) ! Flux reconstructed on each stencil
    real fp    ! Positive flux at cell interfaces
    real fm    ! Negative flux at cell interfaces
    integer i

    umax = abs(u) ! NOTE: In this case, max dflux/drho is just u.

    ! Calculate splitted flux at cell centers by using Lax-Friedrichs splitting.
    ! The splitting is for numerical stability by respecting upwind.
    do i = 1, nx
      fpc(i) = 0.5 * (u + umax) * rho(i)
      fmc(i) = 0.5 * (u - umax) * rho(i)
    end do
    call apply_bc(ns, nx, fpc)
    call apply_bc(ns, nx, fmc)

    do i = 1, nx
      ! Positive flux at cell interfaces
      ! - Calculate flux at interfaces for each stencil.
      fs(1) = dot_product(p1, fpc(i-2:i  ))
      fs(2) = dot_product(p2, fpc(i-1:i+1))
      fs(3) = dot_product(p3, fpc(i  :i+2))
      ! - Calculate smooth indicators for each stencil.
      b(1) = s1 * dot_product(s11, fpc(i-2:i  ))**2 + s2 * dot_product(s12, fpc(i-2:i  ))**2
      b(2) = s1 * dot_product(s21, fpc(i-1:i+1))**2 + s2 * dot_product(s22, fpc(i-1:i+1))**2
      b(3) = s1 * dot_product(s31, fpc(i  :i+2))**2 + s2 * dot_product(s32, fpc(i  :i+2))**2
      ! - Calculate stencil linear combination weights considering smooth indicators.
      w = g / (eps + b)**2
      w = w / sum(w)
      fp = dot_product(w, fs)

      ! Negative flux at cell interfaces
      ! - Calculate flux at interfaces for each stencil.
      fs(1) = dot_product(p1, fmc(i+3:i+1:-1))
      fs(2) = dot_product(p2, fmc(i+2:i  :-1))
      fs(3) = dot_product(p3, fmc(i+1:i-1:-1))
      ! - Calculate smooth indicators for each stencil.
      b(1) = s1 * dot_product(s11, fmc(i+3:i+1:-1))**2 + s2 * dot_product(s12, fmc(i+3:i+1:-1))**2
      b(2) = s1 * dot_product(s21, fmc(i+2:i  :-1))**2 + s2 * dot_product(s22, fmc(i+2:i  :-1))**2
      b(3) = s1 * dot_product(s31, fmc(i+1:i-1:-1))**2 + s2 * dot_product(s32, fmc(i+1:i-1:-1))**2
      ! - Calculate stencil linear combination weights considering smooth indicators.
      w = g / (eps + b)**2
      w = w / sum(w)
      fm = dot_product(w, fs)

      f(i) = fp + fm
    end do
    call apply_bc(1, nx, f)

    dfdx(1:nx) = (f(1:nx) - f(0:nx-1)) / dx

  end subroutine weno

  subroutine rk4()

    ! RK 1st stage
    call weno(rho(:,old), dfdx(:,1))
    rho(1:nx,new) = rho(1:nx,old) - 0.5d0 * dt * dfdx(:,1)
    call apply_bc(ns, nx, rho(:,new))
    ! RK 2nd stage
    call weno(rho(:,new), dfdx(:,2))
    rho(1:nx,new) = rho(1:nx,old) - 0.5d0 * dt * dfdx(:,2)
    call apply_bc(ns, nx, rho(:,new))
    ! RK 3rd stage
    call weno(rho(:,new), dfdx(:,3))
    rho(1:nx,new) = rho(1:nx,old) -         dt * dfdx(:,3)
    call apply_bc(ns, nx, rho(:,new))
    ! RK 4th stage
    call weno(rho(:,new), dfdx(:,4))
    rho(1:nx,new) = rho(1:nx,old) - (dfdx(:,1) + 2.0d0 * dfdx(:,2) + 2.0d0 * dfdx(:,3) + dfdx(:,4)) / 6.0d0
    call apply_bc(ns, nx, rho(:,new))

  end subroutine rk4

  subroutine rk3_tvd()

    ! RK 1st stage
    call weno(rho(:,old), dfdx(:,1))
    rho(1:nx,new) = rho(1:nx,old) - dt * dfdx(1:nx,1)
    call apply_bc(ns, nx, rho(:,new))
    ! RK 2nd stage
    call weno(rho(:,new), dfdx(:,2))
    rho(1:nx,new) = (3.0d0 * rho(1:nx,old) + rho(1:nx,new) - dt * dfdx(:,2)) / 4.0d0
    call apply_bc(ns, nx, rho(:,new))
    ! RK 3rd stage
    call weno(rho(:,new), dfdx(:,3))
    rho(1:nx,new) = (rho(1:nx,old) + 2.0d0 * (rho(1:nx,new) - dt * dfdx(:,3))) / 3.0d0
    call apply_bc(ns, nx, rho(:,new))

  end subroutine rk3_tvd

end program weno_adv_fd_1d_case
