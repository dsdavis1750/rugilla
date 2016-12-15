PRO read_raddat,file,radius,cnts,error

  if(n_params(0) LT 3)then begin
    print,"read_raddat,file,radius,cnts,err"
    print,"Purpose:read a radial plot data file"
    print,"Inputs:"
    print,"  file:name of ds9 radial file"
    print,"Outputs:"
    print,"  radius:vector of radii (pixels)"
    print,"  cnts:vector of cnts/pixel"
    print,"  error:vector of cnts/pixel errors"
    return
   endif

radius = fltarr(35)
cnts = fltarr(35)
error = fltarr(35)
OPENR, lun, file, /GET_LUN
for i = 0, 34 do begin
readf, lun, tmp1, tmp2, tmp3
radius(i) = tmp1
cnts(i) = tmp2
error(i) = tmp3
endfor

FREE_LUN, lun

  end
