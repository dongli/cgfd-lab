! This is a 1D advection example using square initial condition and periodic
! boundary condition for FFSL finite volume scheme.
!
! Li Dong <dongli@lasg.iap.ac.cn>
!
! - 2018-03-22: Initial creation.

program ffsl_adv_1d_case

  use netcdf

  implicit none

  real, allocatable :: x(:)               ! Cell center coordinates
  real, allocatable :: rho(:,:)           ! Tracer density being advected at cell centers
  real, allocatable :: rho_l(:)           ! Tracer density at left cell interfaces
  real, allocatable :: drho(:)            ! Tracer density mismatch at cell centers
  real, allocatable :: rho_6(:)           ! Curvature in PPM at cell centers
  real, allocatable :: flux(:)            ! Flux at cell interfaces
  real dx                                 ! Cell interval
  real :: dt = 1.0                        ! Time step size
  integer :: nx = 100                     ! Cell number
  integer :: nt = 200                     ! Integration time step number
  character(10) :: flux_type = 'ppm'      ! Available flux types: upwind, van_leer, ppm
  character(10) :: limiter_type = 'mono'  ! Available limiter types: none, mono, pd
  real :: u = 0.005                       ! Advection speed
  integer, parameter :: ns = 2            ! Stencil width
  integer i
  integer :: time_step = 0, old = 1, new = 2
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

  allocate(x(nx))
  allocate(rho(1-ns:nx+ns,old:new))
  allocate(rho_l(1-ns:nx+ns))
  allocate(drho(1-ns:nx+ns))
  allocate(rho_6(1-ns:nx+ns))
  allocate(flux(1:nx+1))

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
  print *, time_step, sum(rho(1:nx,old))
  do while (time_step < nt)
    call ffsl(rho(:,old))
    do i = 1, nx
      rho(i,new) = rho(i,old) - (flux(i+1) - flux(i))
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
  deallocate(rho_l)
  deallocate(drho)
  deallocate(rho_6)
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

    x(nx+1) = x(1)

  end subroutine half_boundary_condition

  subroutine ffsl(rho)

    real, intent(in) :: rho(1-ns:nx+ns)

    real cfl, c
    real s1, s2, ds, ds2, ds3
    integer i, j, K, l

    if (flux_type == 'ppm') then
      ! Calculate the subgrid distribution of tracer.
      do i = 1, nx
        call ppm(rho(i-2), rho(i-1), rho(i), rho(i+1), rho(i+2), rho_l(i), drho(i), rho_6(i))
      end do
      call full_boundary_condition(rho_l)
      call full_boundary_condition(drho)
      call full_boundary_condition(rho_6)
    end if

    do i = 1, nx
      flux(i) = 0
      cfl = u * dt / dx
      K = int(cfl)
      c = cfl - K
      ! Calculate integer flux.
      if (cfl > 0) then
        do j = 1, K
          flux(i) = flux(i) + rho(i - j)
        end do
      else if (cfl < 0) then
        do j = 1, -K
          flux(i) = flux(i) - rho(i + j - 1)
        end do
      end if
      l = merge(i - K - 1, i - K, cfl > 0)
      ! Calculate fractional flux.
      select case (flux_type)
      case ('upwind')
        flux(i) = flux(i) + c * rho(l)
      case ('van_leer')
        drho(l) = mismatch(rho(l-1), rho(l), rho(l+1))
        flux(i) = flux(i) + c * (rho(l) + (sign(1.0, cfl) - c) * drho(l) * 0.5)
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
        flux(i) = flux(i) + sign(rho_l(l) * ds + 0.5 * drho(l) * ds2 + rho_6(l) * (0.5 * ds2 - ds3 / 3.0), cfl)
      end select
    end do
    call half_boundary_condition(flux)

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

  subroutine output(rho)

    real, intent(in) :: rho(1-ns:nx+ns)

    character(30) file_name
    integer ncid, time_dimid, time_varid, x_dimid, x_varid, rho_varid, ierr

    write(file_name, "('ffsl.', I3.3, '.', I4.4, '.nc')") nx, time_step

    ierr = NF90_CREATE(file_name, nf90_clobber, ncid)
    ierr = NF90_PUT_ATT(ncid, NF90_GLOBAL, 'scheme', 'FFSL')

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

end program ffsl_adv_1d_case
