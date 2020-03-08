! This is a 1D advection example using square initial condition and periodic
! boundary condition for WENO finite volume scheme.
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
!                                      +|
!                                     i+1/2
!
!   For negative flux:
!
!                       |-----------S1----------|
!                       |       |-----------S2----------|
!                       |       |       |-----------S3----------|
!                    ...|---o---|---o---|---o---|---o---|---o---|...
!                       |  i-1  |   i   |  i+1  |  i+2  |  i+3  |
!                                       |-
!                                     i+1/2
!
! To respect characteristic direction, i.g. upwinding, we need to Calculate
! two versions of density at i+1/2, one is 𝝆+, the other is 𝝆-, then use
! Lax-Friedrichs flux-splitting technique,
!
!   ^           1
!   f(𝝆+,𝝆-) = --- (f(𝝆+) + umax 𝝆+ + f(𝝆-) - umax 𝝆-)
!               2
!
! Li Dong <dongli@lasg.iap.ac.cn>
!
! - 2020-03-04: Initial creation.

program weno_adv_fv_1d_case

  use adv_1d_square_case_mod

  implicit none

  real, allocatable, dimension(:,:) :: rho  ! Tracer density being advected at cell centers
  real, allocatable, dimension(:,:) :: dfdx ! Tendency at cell centers
  real, allocatable, dimension(:  ) :: f    ! Final flux at cell interfaces
  
  real, parameter :: eps = 1.0d-40       ! A small value to avoid divided-by-zero
  integer, parameter :: ns = 3           ! Stencil width
  logical, parameter :: weno_z = .true. ! Switch for using WENO-Z smooth indicator modification


  allocate(rho (1-ns:nx+ns,2))
  allocate(dfdx(1   :nx   ,4))
  allocate(f   (0   :nx+1   ))

  call adv_1d_square_case_init(ns, rho(:,old))
  call output(trim(merge('weno-z', 'weno  ', weno_z)) // '_fv', time_step, ns, nx, x, rho(:,old))

  ! Run integration.
  write(*, *) time_step, sum(rho(1:nx,old))
  do while (time_step < nt)
    call euler()
    ! call rk4()
    ! call rk3_tvd()
    call advance_time()
    call output(trim(merge('weno-z', 'weno  ', weno_z)) // '_fv', time_step, ns, nx, x, rho(:,old))
    write(*, *) time_step, sum(rho(1:nx,old))
  end do

  deallocate(rho)
  deallocate(dfdx)
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

    ! New weight formula from Borges et al. (2008)
    real t5

    real umax    ! Maximum df/drho (e.g. advection velocity in this case)
    real b(3)    ! Smooth indicators
    real w(3)    ! Combination weights
    real rhos(3) ! Density reconstructed on each stencil
    real rhop    ! Positive density at cell interfaces
    real rhom    ! Negative density at cell interfaces
    integer i

    umax = abs(u) ! NOTE: In this case, max dflux/drho is just u.

    do i = 1, nx
      ! Positive flux at cell interfaces
      ! - Calculate flux at interfaces for each stencil.
      rhos(1) = dot_product(p1, rho(i-2:i  ))
      rhos(2) = dot_product(p2, rho(i-1:i+1))
      rhos(3) = dot_product(p3, rho(i  :i+2))
      ! - Calculate smooth indicators for each stencil.
      b(1) = s1 * dot_product(s11, rho(i-2:i  ))**2 + s2 * dot_product(s12, rho(i-2:i  ))**2
      b(2) = s1 * dot_product(s21, rho(i-1:i+1))**2 + s2 * dot_product(s22, rho(i-1:i+1))**2
      b(3) = s1 * dot_product(s31, rho(i  :i+2))**2 + s2 * dot_product(s32, rho(i  :i+2))**2
      ! - WENO-Z
      if (weno_z) then
        t5 = abs(b(1) - b(3))
        b = (b + eps) / (b + t5 + eps)
      end if
      ! - Calculate stencil linear combination weights considering smooth indicators.
      w = g / (eps + b)**2
      w = w / sum(w)
      rhop = dot_product(w, rhos)

      ! Negative flux at cell interfaces
      ! - Calculate flux at interfaces for each stencil.
      rhos(1) = dot_product(p1, rho(i+3:i+1:-1))
      rhos(2) = dot_product(p2, rho(i+2:i  :-1))
      rhos(3) = dot_product(p3, rho(i+1:i-1:-1))
      ! - Calculate smooth indicators for each stencil.
      b(1) = s1 * dot_product(s11, rho(i+3:i+1:-1))**2 + s2 * dot_product(s12, rho(i+3:i+1:-1))**2
      b(2) = s1 * dot_product(s21, rho(i+2:i  :-1))**2 + s2 * dot_product(s22, rho(i+2:i  :-1))**2
      b(3) = s1 * dot_product(s31, rho(i+1:i-1:-1))**2 + s2 * dot_product(s32, rho(i+1:i-1:-1))**2
      ! - WENO-Z
      if (weno_z) then
        t5 = abs(b(1) - b(3))
        b = (b + eps) / (b + t5 + eps)
      end if
      ! - Calculate stencil linear combination weights considering smooth indicators.
      w = g / (eps + b)**2
      w = w / sum(w)
      rhom = dot_product(w, rhos)

      f(i) = 0.5 * ((u + umax) * rhop + (u - umax) * rhom)
    end do
    call apply_bc(1, nx, f)

    dfdx(1:nx) = (f(1:nx) - f(0:nx-1)) / dx

  end subroutine weno

  subroutine euler()

    call weno(rho(:,old), dfdx(:,1))
    rho(1:nx,new) = rho(1:nx,old) - dt * dfdx(:,1)
    call apply_bc(ns, nx, rho(:,new))

  end subroutine euler

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

end program weno_adv_fv_1d_case
