FUNCTION xmm_off_axis_ang,detx,dety,inst

  ;program to determine the distance from the optical axis in arcmin
  ;given the detector coordinates

  case inst of
  "EMOS1":begin
    radi=sqrt((detx-100)^2+(dety+200)^2)
  end
  "EMOS2":begin
    radi=sqrt((detx-500)^2+(dety+1250)^2)
  end
  "EPN":begin
    radi=sqrt((detx-1240)^2+(dety-400)^2)
  end
  endcase

  radi=radi*0.05/60. ;answer returned in arcmin
  return,radi
  end


