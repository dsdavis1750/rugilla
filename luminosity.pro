
;solar luminosity=3.9*10^33 erg/s

z_soft=[0.055, 0.058, 0.039, 0.06, 0.057, 0.06, 0.06, 0.06, 0.06, $
	0.065, 0.06, 0.063, 0.06]
B_soft=[15.94, 17.63, 15.28, 21.28, 15.82, 18.22, 20.09, 20.42, 18.38, 16.26, $
	20.48, 18.08, 20.32]
optical_soft=2/5.*(0.48-B_soft+5*alog10(4.63*10^9.*((1+z_soft)^2.-1)/((1+z_soft)^2.+1)))

xray_soft=[0.55, 2.9, 0.56, 0.94, 4.3, 0.56, 8.1, 17, 0.96, 1.7, 0.68, 0.23]
	;in units of 10^41 erg/s
rate_soft=[1.61, 7.3, 3.19, 2.48, 10.7, 1.3, 20.0, 42.4, 2.03, 3.93, 1.54, 0.568]
error_soft=xray_soft/sqrt(19.7632*rate_soft)

optical_soft7=10.8	;Source #7 info
xray_soft7=3.3
error_soft7=0.25


z_hard=[0.058, 0.06, 0.06, 0.06, 0.06, 0.06, 0.06, 0.063]
B_hard=[17.63, 21.28, 18.22, 20.09, 20.42, 18.38, 20.48, 18.08]
optical=2/5.*(0.48-B_hard+5*alog10(4.63*10^9.*((1+z_hard)^2.-1)/((1+z_hard)^2.+1)))
optical_hard=[9.5, 9.7, 9.1, 9.8, 10.5]
xray_hard=[7.5, 1.3, 5.6, 1.2, 13, 31, 2.0, 2.1]
	;in units of 10^41 erg/s
error_hard=[0.81, 0.0, 0.71, 0.0, 1.1, 1.6, 0.0, 0.46]

optical_hard7=10.8	;Source #7 info
xray_hard7=0.79
error_hard7=0.25

;open_plot, 'soft.ps'

;plot, optical_soft, 41+alog10(xray_soft), charsize=1.5, $;
;psym=7, yrange=[38, 44], xtitle='log L!LB!N/L!LSun!UB!N', $
;ytitle='log L!Lx!N(soft) (erg/s)', xrange=[8.5, 11]
;errplot, optical_soft, 41+alog10(xray_soft)+alog10(1-error_soft/xray_soft), $
;41+alog10(xray_soft)+alog10(1+error_soft/xray_soft)
;plots, [8.5,11],[38.1,40.6], linestyle=2
;plots, optical_soft7, 41+alog10(xray_soft7), psym=7
;errplot, optical_soft7, $
;41+alog10(xray_soft7)+alog10(1-error_soft7/xray_soft7), $
;41+alog10(xray_soft7)+alog10(1+error_soft7/xray_soft7)

;plot, 10-alog10(1.95)+alog10(optical_hard), 41+alog10(xray_hard), $
plot, optical_hard, 41+alog10(xray_hard), charsize=1.7, $
xrange=[8.5,11], yrange=[38,44], psym=7, xtitle='!5log L!LB!N/L!LSun!UB!N', $
ytitle='L!Lx!N(hard) (10!U41!Nerg/s)'
errplot, optical_hard, 41+alog10(xray_hard)+alog10(1-error_hard/xray_hard), $
41+alog10(xray_hard)+alog10(1+error_hard/xray_hard)
plots, [8.5,11],[38.1,40.6], linestyle=2
plots, optical_hard7, 41+alog10(xray_hard7), psym=2
errplot, optical_hard7, $
41+alog10(xray_hard7)+alog10(1-error_hard7/xray_hard7), $
41+alog10(xray_hard7)+alog10(1+error_hard7/xray_hard7)

;close_plot

;write_gif, 'my_soft.gif', tvrd()
end
