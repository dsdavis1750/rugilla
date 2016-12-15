pro areas_p,dummy
;
;
tbread,'test_2',h2,tab2
areas=tbget(h2,tab2,1)
!xmin=10.0 & !xmax=729
!ymin= 0.01 & !ymax=1000
!xtitle='channel #'
!ytitle='effective area (cm**2)'
plot_oo,areas(0,*)
oplot,areas(1,*)
oplot,areas(2,*)
oplot,areas(3,*)
oplot,areas(4,*)
oplot,areas(5,*)
oplot,areas(6,*)
oplot,areas(7,*)
oplot,areas(8,*)
oplot,areas(9,*)
oplot,areas(10,*)
oplot,areas(11,*)
oplot,areas(12,*)
oplot,areas(13,*)
return
end

