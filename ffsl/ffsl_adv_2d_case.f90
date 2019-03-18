! This is a 2D advection example using slotted cylinder initial condition and
! double periodic boundary condition for FFSL finite volume scheme.
!
! Li Dong <dongli@lasg.iap.ac.cn>
!
! - 2018-03-24: Initial creation.
! - 2018-03-25: Remove flux splitting (i.e. integer and fractional flux parts),
!               because in 2D this is not very helpful to increase stability,
!               we are still limited by CFL condition.
! - 2019-03-18: Improve readability by add inner and outer operators comments
!               and add flux_x and flux_y arguments to subroutine ffsl.

program ffsl_adv_2d_case

  use netcdf

  implicit none

  real, allocatable :: x(:)               ! Cell center coordinates
  real, allocatable :: y(:)               ! Cell center coordinates
  real, allocatable :: u(:,:)             ! Velocity component along x axis
  real, allocatable :: v(:,:)             ! Velocity component along y axis
  real, allocatable :: rho(:,:,:)         ! Tracer density being advected at cell centers
  real, allocatable :: rho_x(:,:)         ! Tracer density due to advective operator along x axis
  real, allocatable :: rho_y(:,:)         ! Tracer density due to advective operator along y axis
  real, allocatable :: rho_xl(:,:)        ! Tracer density at left cell interfaces along x axis
  real, allocatable :: rho_yl(:,:)        ! Tracer density at left cell interfaces along y axis
  real, allocatable :: drho_x(:,:)        ! Tracer density mismatch at cell centers along x axis
  real, allocatable :: drho_y(:,:)        ! Tracer density mismatch at cell centers along y axis
  real, allocatable :: rho_x6(:,:)        ! PPM polynomial coefficient at cell centers along x axis
  real, allocatable :: rho_y6(:,:)        ! PPM polynomial coefficient at cell centers along y axis
  real, allocatable :: flux_x(:,:)        ! Flux at cell interfaces along x axis
  real, allocatable :: flux_y(:,:)        ! Flux at cell interfaces along y axis
  real, allocatable :: div_x(:,:)         ! Divergence component along x axis
  real, allocatable :: div_y(:,:)         ! Divergence component along y axis
  real dx                                 ! Cell interval along x axis
  real dy                                 ! Cell interval along y axis
  real :: dt = 0.25                       ! Time step size
  integer :: nx = 100                     ! Cell number along x axis
  integer :: ny = 100                     ! Cell number along y axis
  integer :: nt = 836                     ! Integration time step number
  character(10) :: flux_type = 'ppm'      ! Available flux types: upwind, van_leer, ppm
  character(10) :: limiter_type = 'mono'  ! Available limiter types: none, mono, pd
  real, parameter :: omega = 0.03         ! Rotation angular speed
  real, parameter :: x0 = 0.25            ! Initial coordinate
  real, parameter :: y0 = 0.5             ! Initial coordinate
  integer, parameter :: ns = 2            ! Stencil width
  integer i, j
  integer :: time_step = 0, old = 1, new = 2
  character(256) namelist_path
  logical is_exist

  namelist /params/ nx, ny, nt, dt, flux_type, limiter_type

  call get_command_argument(1, namelist_path)
  inquire(file=namelist_path, exist=is_exist)
  if (is_exist) then
    open(10, file=namelist_path)
    read(10, nml=params)
    close(10)
  end if

  allocate(x(nx))
  allocate(y(ny))
  allocate(u(1:nx+1,ny))
  allocate(v(nx,1:ny+1))
  allocate(rho(1-ns:nx+ns,1-ns:ny+ns,old:new))
  allocate(rho_x(1-ns:nx+ns,1-ns:ny+ns))
  allocate(rho_y(1-ns:nx+ns,1-ns:ny+ns))
  allocate(rho_xl(1-ns:nx+ns,1-ns:ny+ns))
  allocate(rho_yl(1-ns:nx+ns,1-ns:ny+ns))
  allocate(drho_x(1-ns:nx+ns,1-ns:ny+ns))
  allocate(drho_y(1-ns:nx+ns,1-ns:ny+ns))
  allocate(rho_x6(1-ns:nx+ns,1-ns:ny+ns))
  allocate(rho_y6(1-ns:nx+ns,1-ns:ny+ns))
  allocate(flux_x(1:nx+1,ny))
  allocate(flux_y(nx,1:ny+1))
  allocate(div_x(1-ns:nx+ns,1-ns:ny+ns))
  allocate(div_y(1-ns:nx+ns,1-ns:ny+ns))

  ! Set mesh grid coordinates.
  dx = 1.0d0 / nx
  do i = 1, nx
    x(i) = (i - 1) * dx
  end do
  dy = 1.0d0 / ny
  do j = 1, ny
    y(j) = (j - 1) * dy
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
  call half_x_boundary_condition(u)
  call half_y_boundary_condition(v)
  call full_boundary_condition(rho(:,:,old))
  call output(rho(:,:,old))

  ! Run integration.
  print *, time_step, sum(rho(1:nx,1:ny,old))
  do while (time_step < nt)
    ! --------------------------------------------------------------------------
    ! Run inner operators.
    call ffsl(rho(:,:,old), rho(:,:,old), flux_x, flux_y)
    call divergence()
    ! Subtract divergence terms from flux to form advective operators.
    do j = 1, ny
      do i = 1, nx + 1
        flux_x(i,j) = flux_x(i,j) - 0.5 * (div_x(i,j) * rho(i,j,old) + div_x(i-1,j) * rho(i-1,j,old))
      end do
    end do
    do j = 1, ny + 1
      do i = 1, nx
        flux_y(i,j) = flux_y(i,j) - 0.5 * (div_y(i,j) * rho(i,j,old) + div_y(i,j-1) * rho(i,j-1,old))
      end do
    end do
    ! Calculate intermediate tracer density due to advective operators.
    do j = 1, ny
      do i = 1, nx
        rho_x(i,j) = rho(i,j,old) - 0.5 * (flux_x(i+1,j) - flux_x(i,j))
        rho_y(i,j) = rho(i,j,old) - 0.5 * (flux_y(i,j+1) - flux_y(i,j))
      end do
    end do
    call full_boundary_condition(rho_x)
    call full_boundary_condition(rho_y)
    ! --------------------------------------------------------------------------
    ! Run outer operators.
    call ffsl(rho_y, rho_x, flux_x, flux_y)
    do j = 1, ny
      do i = 1, nx
        rho(i,j,new) = rho(i,j,old) - (flux_x(i+1,j) - flux_x(i,j)) - (flux_y(i,j+1) - flux_y(i,j))
      end do
    end do
    call full_boundary_condition(rho(:,:,new))
    ! Change time indices.
    i = old; old = new; new = i
    time_step = time_step + 1
    call output(rho(:,:,old))
    print *, time_step, sum(rho(1:nx,1:ny,old))
  end do

  deallocate(x)
  deallocate(rho)
  deallocate(rho_x)
  deallocate(rho_y)
  deallocate(rho_xl)
  deallocate(rho_yl)
  deallocate(drho_x)
  deallocate(drho_y)
  deallocate(rho_x6)
  deallocate(rho_y6)
  deallocate(flux_x)
  deallocate(flux_y)
  deallocate(div_x)
  deallocate(div_y)

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

  subroutine half_x_boundary_condition(array)

    real, intent(inout) :: array(1:nx+1,ny)

    array(nx+1,:) = array(1,:)

  end subroutine half_x_boundary_condition

  subroutine half_y_boundary_condition(array)

    real, intent(inout) :: array(nx,1:ny+1)

    array(:,ny+1) = array(:,1)

  end subroutine half_y_boundary_condition

  subroutine ffsl(rho_x, rho_y, flux_x, flux_y)

    real, intent(in) :: rho_x(1-ns:nx+ns,1-ns:ny+ns)
    real, intent(in) :: rho_y(1-ns:nx+ns,1-ns:ny+ns)
    real, intent(out) :: flux_x(1:nx+1,ny)
    real, intent(out) :: flux_y(nx,1:ny+1)

    real c, s1, s2, ds, ds2, ds3, max_cfl
    integer i, j, k

    if (flux_type == 'ppm') then
      ! Calculate the subgrid distribution of tracer.
      do j = 1, ny
        do i = 1, nx
          call ppm(rho_x(i-2,j), rho_x(i-1,j), rho_x(i,j), rho_x(i+1,j), rho_x(i+2,j), rho_xl(i,j), drho_x(i,j), rho_x6(i,j))
          call ppm(rho_y(i,j-2), rho_y(i,j-1), rho_y(i,j), rho_y(i,j+1), rho_y(i,j+2), rho_yl(i,j), drho_y(i,j), rho_y6(i,j))
        end do
      end do
      call full_boundary_condition(rho_xl)
      call full_boundary_condition(rho_yl)
      call full_boundary_condition(drho_x)
      call full_boundary_condition(drho_y)
      call full_boundary_condition(rho_x6)
      call full_boundary_condition(rho_y6)
    end if

    max_cfl = -1
    ! Along x axis
    do j = 1, ny
      do i = 1, nx
        c = u(i,j) * dt / dx
        max_cfl = max(c, max_cfl)
        k = merge(i - 1, i, c > 0)
        select case (flux_type)
        case ('upwind')
          flux_x(i,j) = c * rho_x(k,j)
        case ('van_leer')
          drho_x(k,j) = mismatch(rho_x(k-1,j), rho_x(k,j), rho_x(k+1,j))
          flux_x(i,j) = c * (rho_x(k,j) + (sign(1.0, c) - c) * drho_x(k,j) * 0.5d0)
        case ('ppm')
          if (c >= 0) then
            s1 = 1 - abs(c)
            s2 = 1
          else if (c < 0) then
            s1 = 0
            s2 = abs(c)
          end if
          ds = s2 - s1
          ds2 = s2**2 - s1**2
          ds3 = s2**3 - s1**3
          flux_x(i,j) = sign(rho_xl(k,j) * ds + 0.5 * drho_x(k,j) * ds2 + rho_x6(k,j) * (0.5 * ds2 - ds3 / 3.0), c)
        end select
      end do
    end do
    call half_x_boundary_condition(flux_x)
    ! Along y axis
    do j = 1, ny
      do i = 1, nx
        c = v(i,j) * dt / dy
        max_cfl = max(c, max_cfl)
        k = merge(j - 1, j, c > 0)
        select case (flux_type)
        case ('upwind')
          flux_y(i,j) = c * rho_y(i,k)
        case ('van_leer')
          drho_y(i,k) = mismatch(rho_y(i,k-1), rho_y(i,k), rho_y(i,k+1))
          flux_y(i,j) = c * (rho_y(i,k) + (sign(1.0, c) - c) * drho_y(i,k) * 0.5d0)
        case ('ppm')
          if (c >= 0) then
            s1 = 1 - abs(c)
            s2 = 1
          else if (c < 0) then
            s1 = 0
            s2 = abs(c)
          end if
          ds = s2 - s1
          ds2 = s2**2 - s1**2
          ds3 = s2**3 - s1**3
          flux_y(i,j) = sign(rho_yl(i,k) * ds + 0.5 * drho_y(i,k) * ds2 + rho_y6(i,k) * (0.5 * ds2 - ds3 / 3.0), c)
        end select
      end do
    end do
    call half_y_boundary_condition(flux_y)
    print *, max_cfl

  end subroutine ffsl

  subroutine divergence()

    integer i, j

    do j = 1, ny
      do i = 1, nx
        div_x(i,j) = dt / dx * (u(i+1,j) - u(i,j))
        div_y(i,j) = dt / dy * (v(i,j+1) - v(i,j))
      end do
    end do
    call full_boundary_condition(div_x)
    call full_boundary_condition(div_y)

  end subroutine divergence

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

    real, intent(in) :: rho(1-ns:nx+ns,1-ns:ny+ns)

    character(30) file_name
    integer file_id, ierr
    integer time_dim_id, time_var_id
    integer x_dim_id, x_var_id
    integer y_dim_id, y_var_id
    integer rho_var_id

    write(file_name, "('ffsl.', I3.3, '.nc')") time_step

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

    ierr = nf90_put_var(file_id, x_var_id, x)
    if (ierr /= nf90_noerr) then
      write(*, *) '[Error]: Failed to write x variable!'
      stop 1
    end if

    ierr = nf90_put_var(file_id, y_var_id, y)
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

end program ffsl_adv_2d_case
