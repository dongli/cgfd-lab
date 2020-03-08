! This is a 1D advection example using square initial condition and periodic
! boundary condition for FFSL finite volume scheme.
!
! Li Dong <dongli@lasg.iap.ac.cn>
!
! - 2018-03-22: Initial creation.

program ffsl_adv_1d_case

  use adv_1d_square_case_mod

  implicit none

  real, allocatable :: rho(:,:)           ! Tracer density being advected at cell centers
  real, allocatable :: rhol(:)            ! Tracer density at left cell interfaces
  real, allocatable :: drho(:)            ! Tracer density mismatch at cell centers
  real, allocatable :: rho6(:)            ! Curvature in PPM at cell centers
  real, allocatable :: f(:)               ! Flux at cell interfaces
  character(10) :: flux_type = 'ppm'      ! Available flux types: upwind, van_leer, ppm
  character(10) :: limiter_type = 'mono'  ! Available limiter types: none, mono, pd
  integer, parameter :: ns = 2            ! Stencil width
  integer i
  character(256) namelist_path
  logical is_exist

  namelist /params/ nx, nt, dt, flux_type, limiter_type, u

  call get_command_argument(1, namelist_path)
  inquire(file=namelist_path, exist=is_exist)
  if (is_exist) then
    open(10, file=namelist_path)
    read(10, nml=params)
    close(10)
  end if

  allocate(rho (1-ns:nx+ns,2))
  allocate(rhol(1-ns:nx+ns))
  allocate(drho(1-ns:nx+ns))
  allocate(rho6(1-ns:nx+ns))
  allocate(f   (   0:nx+1 ))

  call adv_1d_square_case_init(ns, rho(:,old))
  call output('ffsl', time_step, ns, nx, x, rho(:,old))

  ! Run integration.
  print *, time_step, sum(rho(1:nx,old))
  do while (time_step < nt)
    call ffsl(rho(:,old))
    do i = 1, nx
      rho(i,new) = rho(i,old) - (f(i+1) - f(i))
    end do
    call apply_bc(ns, nx, rho(:,new))
    call advance_time()
    call output('ffsl', time_step, ns, nx, x, rho(:,old))
    print *, time_step, sum(rho(1:nx,old))
  end do

  deallocate(rho)
  deallocate(rhol)
  deallocate(drho)
  deallocate(rho6)
  deallocate(f)

  call adv_1d_square_case_final()

contains

  subroutine ffsl(rho)

    real, intent(in) :: rho(1-ns:nx+ns)

    real cfl, c
    real s1, s2, ds, ds2, ds3
    integer i, j, K, l

    if (flux_type == 'ppm') then
      ! Calculate the subgrid distribution of tracer.
      do i = 1, nx
        call ppm(rho(i-2), rho(i-1), rho(i), rho(i+1), rho(i+2), rhol(i), drho(i), rho6(i))
      end do
      call apply_bc(ns, nx, rhol)
      call apply_bc(ns, nx, drho)
      call apply_bc(ns, nx, rho6)
    end if

    do i = 1, nx
      f(i) = 0
      cfl = u * dt / dx
      K = int(cfl)
      c = cfl - K
      ! Calculate integer flux.
      if (cfl > 0) then
        do j = 1, K
          f(i) = f(i) + rho(i - j)
        end do
      else if (cfl < 0) then
        do j = 1, -K
          f(i) = f(i) - rho(i + j - 1)
        end do
      end if
      l = merge(i - K - 1, i - K, cfl > 0)
      ! Calculate fractional flux.
      select case (flux_type)
      case ('upwind')
        f(i) = f(i) + c * rho(l)
      case ('van_leer')
        drho(l) = mismatch(rho(l-1), rho(l), rho(l+1))
        f(i) = f(i) + c * (rho(l) + (sign(1.0, cfl) - c) * drho(l) * 0.5)
      case ('ppm')
        if (cfl >= 0) then
          s1 = 1 - abs(c)
          s2 = 1
        else if (cfl < 0) then
          s1 = 0
          s2 = abs(c)
        end if
        ds = s2 - s1
        ds2 = s2**2 - s1**2
        ds3 = s2**3 - s1**3
        f(i) = f(i) + sign(rhol(l) * ds + 0.5 * drho(l) * ds2 + rho6(l) * (0.5 * ds2 - ds3 / 3.0), cfl)
      end select
    end do
    call apply_bc(1, nx, f)

  end subroutine ffsl

  subroutine ppm(fm2, fm1, f, fp1, fp2, fl, df, f6)

    real, intent(in) :: fm2
    real, intent(in) :: fm1
    real, intent(in) :: f
    real, intent(in) :: fp1
    real, intent(in) :: fp2
    real, intent(out) :: fl
    real, intent(out) :: df
    real, intent(out) :: f6

    real dfl, dfr, fr

    ! Calculate values at left and right cell interfaces.
    dfl = mismatch(fm2, fm1, f  )
    df  = mismatch(fm1, f,   fp1)
    dfr = mismatch(f,   fp1, fp2)
    ! Why (B2) in Lin (2004) divide (dfl - df) and (df - dfr) by 3?
    fl = 0.5 * (fm1 + f) + (dfl - df) / 6.0
    fr = 0.5 * (fp1 + f) + (df - dfr) / 6.0
    ! Why (B3) and (B4) in Lin (2004) multiply df by 2?
    fl = f - sign(min(abs(df), abs(fl - f)), df)
    fr = f + sign(min(abs(df), abs(fr - f)), df)
    f6 = 6 * f - 3 * (fl + fr)
    df = fr - fl

  end subroutine ppm

  real function mismatch(fm1, f, fp1)

    real, intent(in) :: fm1
    real, intent(in) :: f
    real, intent(in) :: fp1

    real df, df_min, df_max

    df = (fp1 - fm1) * 0.5
    select case (limiter_type)
    case ('none')
      mismatch = df
    case ('mono')
      df_min = 2 * (f - min(fm1, f, fp1))
      df_max = 2 * (max(fm1, f, fp1) - f)
      mismatch = sign(min(abs(df), df_min, df_max), df)

      ! The following codes are (1.8) from Collela and Woodward (1984). It should be equivalent with the above.
      ! if ((fp1 - f) * (f - fm1) > 0) then
      !   mismatch = sign(min(abs(df), abs(f - fm1), abs(f - fp1)), df)
      ! else
      !   mismatch = 0.0
      ! end if
    case ('pd')
      mismatch = sign(min(abs(df), 2 * f), df)
    end select

  end function mismatch

end program ffsl_adv_1d_case
