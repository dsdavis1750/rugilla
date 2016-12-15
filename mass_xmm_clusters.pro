pro mass_xmm,dummy
;
; physical constants
   G = 6.67d-8
   k = 1.38d-16            ; erg/K
   Mpc = 3.08d24           ; cm
   M_sun = 1.989d33        ; grams
   proton_mass = 1.67d-24  ; grams
   erg_per_ev = 1.602d-12
   mu = 0.62               ; mean molecular weight for fully ionized gas
;
; Observational parameters
;
rcore = 0.4694    ; Mpc
rmax = 3.0e-1     ; Mpc
temp = 5.12       ; keV
beta = 0.730      ; slope of surface brightness
n0 = 2.03e-3      ; central electron number density
H     = 70   ; km/s/Mpc
scale = 3.26 ; arcsec/kpc
h_scale = H/100.

;
; radial vector
x = findgen(3000) + 1.
x = (x / 30.)      ; set range from 0 - 30 rcore
;
;  Mvir from Cowie et al. 1986?
gamma = 1.
delta = 3./2 * beta /( 1. - 0.2*(gamma - 1.))
phi = delta * ( gamma - 1.)
a = rcore
t10 = temp/10.
mvir = 1.92e14 * delta * gamma * t10 * (a/0.25)*(x^3)*(1+x^2)^(-(1.+phi))
mgas = 4.5e13*(a/0.4)^3*n0/.002 * (x/2)*sqrt(x*x + 1.) - 0.5*(alog10(x + sqrt(x*x+1.)))
;
;
dtdx=0.
r = x * rcore * Mpc ; radius in cm
kT = (temp * 1000.0) * erg_per_ev
sigma_sq = kT/(mu*proton_mass)
x_factor = (x^3) / (1 + x^2)
m_tot = (3*beta*sigma_sq*x_factor*rcore*Mpc)/G - (x^2 *k* sigma_sq/(G*kT))*dtdx 
m_tot = m_tot / M_sun ; convert to solar masses
;
; electron distribution
n_elec = n0 * ( (1. + x^2))^(-1.5 * beta)
;
entropy = (temp / n_elec^(2./3.))
;
;
x_plot = x*rcore*1000. ; plot in kpc
!xtitle = 'Radius (kpc)'
;
r_vir = 1.9*sqrt(temp/10)/h_scale ; Thomas et al. 2001 MNRAS 324 450
print,'Virial radius: ',r_vir, " Mpc"
x_plot = x*rcore*1000./r_vir ; plot in units of R_vir
!xtitle = 'Radius (r/r!dvir!n)'
;
;
arr_end = n_elements(m_tot)-1
print,'Total Mass:',m_tot(arr_end),' at radius: ',r(arr_end)/Mpc
print,'Total Mass:',m_tot(63),' at radius: ',r(63)/Mpc
;
x_plot = x ; plot in units of core radii
;!xtitle = 'Radius (r/r!dcore!n)'
;
; plotting stuff
;
;
!p.multi=[0,2,2]
!x.style=3
!y.style=3
x_max = max(x_plot)/4.
x_max = x_plot(63)
x_min = min(x_plot)
;
;
!ytitle = 'Mass (Solar Masses)'
!p.title = 'Mass profile - Abell 222'
!p.linestyle = 0
;plot_oo,x*rcore*1000.,m_tot,yrange=[1.e7,1.e16],xrange=[1.0,300.]
plot_oo,x_plot,m_tot,yrange=[1.e7,1.e16],xrange=[x_min,x_max]
!p.linestyle = 3
oplot,x_plot,mgas
;
; Entropy plot
!p.linestyle = 0
!p.title = 'Entropy profile - Abell 222'
!ytitle='Entropy (keV cm!u2!d)'
;plot_oi,x_plot,entropy,yrange=[30.,250.],xrange=[1.0,300.]
plot_oo,x_plot,entropy,yrange=[1.,10000.],xrange=[x_min,x_max]
;plot,x_plot,entropy,yrange=[10.,150.],xrange=[1.0,3.00]
;
; Density plot
!p.linestyle = 0
!p.title ='Density profile - Abell 222'
!ytitle='n!de!n (cm!u-3!d)'
plot_oo,x_plot,n_elec,yrange=[1.e-6,4.e-2],xrange=[x_min,x_max]
;
; temperature plot
!p.title ='Temperature profile - Abell 222'
!ytitle='T(keV)'
plot_oi,x_plot,fltarr(n_elements(n_elec))+temp,yrange=[5.0,6.0],xrange=[x_min,x_max]
!p.multi=0
;print,entropy
;pprint,kT/(1000.*erg_per_ev),temp,x(999)
indx = where((x gt 0.0999) and (x lt 0.109))
;print,x(indx),entropy(indx)
;stop
return
end
