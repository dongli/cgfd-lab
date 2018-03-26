load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

  available_schemes = (/ \
    (/"ftcs",           "Forward-time center-space scheme"/), \
    (/"upwind",         "Upwind scheme"/), \
    (/"beam_warming",   "Beam-Warming scheme"/), \
    (/"lax_wendroff",   "Lax-Wendroff scheme"/), \
    (/"fromm",          "Fromm scheme"/), \
    (/"leap_frog",      "Leap-frog scheme"/), \
    (/"crank_nicolson", "Crank-Nicolson scheme"/), \
    (/"tspas",          "TSPAS scheme"/), \
    (/"mpdata",         "MPDATA scheme"/), \
    (/"ffsl",           "FFSL scheme"/), \
    (/"weno",           "WENO scheme"/) \
  /)

  if (isdefined("scheme")) then
    dims = dimsizes(available_schemes)
    do i = 0, dims(0) - 1
      if (available_schemes(i,0) .eq. scheme) then
        scheme_index = i
        break
      end if
    end do
  else
    print("[Error]: You need select a scheme to plot!")
    exit
  end if

  errors = new(4, float)

  nx = (/100, 200, 400, 800/)
  nt = (/200, 400, 800, 1600/)
  do i = 0, dimsizes(nx) - 1
    f0 = addfile(available_schemes(scheme_index, 0) + "." + nx(i) + ".0000.nc", "r")
    f1 = addfile(available_schemes(scheme_index, 0) + "." + nx(i) + "." + sprinti("%4.4d", nt(i)) + ".nc", "r")

    dx = f0->x(1) - f0->x(0)
    rho0 := f0->rho(0,:)
    rho1 := f1->rho(0,:)
    diff := rho1 - rho0
    errors(i) = sqrt(sum(diff * diff * dx) / sum(rho0 * rho0 * dx))
  end do
  print(errors)

  wks = gsn_open_wks("pdf", available_schemes(scheme_index, 0) + ".errors")

  res = True
  res@tiMainString = available_schemes(scheme_index, 1)
  res@gsnDraw = False
  res@gsnFrame = False
  res@tiXAxisString = "Grid Number"
  res@tiYAxisString = "L2 Error"
  res@xyExplicitLegendLabels = "" + nx
  res@xyLineThicknessF = 5.0
  res@xyYStyle = "Log"
  res@xyComputeYMin = True
  res@xyXStyle = "Log"
  res@pmTickMarkDisplayMode = True
  res@tmXBMode = "Explicit"
  res@tmXBValues = nx
  res@tmXBLabels = "" + nx

  plot = gsn_csm_xy(wks, nx, errors, res)

  ref1 = (1.0 / nx)
  ref2 = 20 * (1.0 / nx)^2
  ref3 = 2 * (1.0 / nx)^3

  res@xyLineColor = "red"
  plot_ref1 = gsn_csm_xy(wks, nx, ref1, res)
  res@xyLineColor = "blue"
  plot_ref2 = gsn_csm_xy(wks, nx, ref2, res)
  res@xyLineColor = "green"
  plot_ref3 = gsn_csm_xy(wks, nx, ref3, res)

  overlay(plot, plot_ref1)
  overlay(plot, plot_ref2)
  overlay(plot, plot_ref3)

  draw(plot)
  frame(wks)

end
