PRO dops,ori,col,xsiz,ysiz

  set_plot,'PS'
  if(ori EQ 'p' or ori EQ 'P')then begin
    device,/portrait,xoffset=3,yoffset=2,xsize=xsiz,ysize=ysiz
  endif else begin
    device,/landscape,xoffset=3,yoffset=27,xsize=xsiz,ysize=ysiz
  endelse

  if(col)then begin
    device,/color,bits_per_pixel=8
  endif

  end
