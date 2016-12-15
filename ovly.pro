PRO ovly,a,x,y,lvl
SET_PLOT,1
save_psym = !psym
!NOERAS=-1 & !psym=0
S=SIZE(A)
if s(0) ne 2 then begin
	print,'image not 2d'
	return 
	endif
set_screen,x,x+s(1)-1,y,y+s(2)-1
set_xy,0,s(1)-1,0,s(2)-1
!type=12
contour,a,lvl
!type=0 & !noeras=0
!ymax=0 &set_plot,0
!psym = save_psym
return 
end
