pro mass_n4325,dummy
;
; physical constants
   G = 6.67d-8
   k = 1.38d-16            ; erg/K
   Mpc = 3.08d24           ; cm
   M_sun = 1.989d33        ; grams
   proton_mass = 1.67d-24  ; grams
   erg_per_ev = 1.602d-12
;
; Observational parameters
;
rcore = 1.53e-2  ; Mpc
rmax = 3.0e-1     ; Mpc
temp = 0.895      ; keV
mu = 0.62         ; mean molecular weight for fully ionized gas
beta = 0.560      ; slope of surface brightness
n0 = 2.09e-2      ; central electron number density
;
; radial vector
x = findgen(3000) + 1.
x = (x / 30.)      ; set range from 0 - 30 rcore
;x = x /rcore        ; convert x to units of rcore
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
x_plot = x*rcore*1000./1080 ; plot in units of R_vir
!xtitle = 'Radius (r/r!dvir!n)'
;
;x_plot = x ; plot in units of core radii
;!xtitle = 'Radius (r/r!dcore!n)'
;
; plotting stuff
;
;
!p.multi=[0,2,2]
!x.style=3
!y.style=3
x_max = max(x_plot)/4.
x_min = min(x_plot)
;
;
!ytitle = 'Mass (Solar Masses)'
!p.title = 'Mass profile - NGC 4325 Group'
!p.linestyle = 0
;plot_oo,x*rcore*1000.,m_tot,yrange=[1.e7,1.e14],xrange=[1.0,300.]
plot_oo,x_plot,m_tot,yrange=[1.e7,1.e14],xrange=[x_min,x_max]
!p.linestyle = 8
oplot,x_plot,mgas
;
; Entropy plot
!p.linestyle = 0
!p.title = 'Entropy profile - NGC 4325 Group'
!ytitle='Entropy (keV cm!u2!d)'
;plot_oi,x_plot,entropy,yrange=[30.,250.],xrange=[1.0,300.]
plot,x_plot,entropy,yrange=[10.,450.],xrange=[x_min,x_max]
;plot,x_plot,entropy,yrange=[10.,150.],xrange=[1.0,3.00]
;
; Density plot
!p.linestyle = 0
!p.title ='Density profile - NGC 4325 Group'
!ytitle='n!de!n (cm!u-3!d)'
plot_oo,x_plot,n_elec,yrange=[1.e-4,4.e-2],xrange=[x_min,x_max]
;
; temperature plot
!p.title ='Temperature profile - NGC 4325 Group'
!ytitle='T(keV)'
plot,x_plot,fltarr(n_elements(n_elec))+0.95,yrange=[0.9,1.0],xrange=[x_min,x_max]
!p.multi=0
print,entropy
print,kT/(1000.*erg_per_ev),temp,x(999)
indx = where((x gt 0.0999) and (x lt 0.109))
print,x(indx),entropy(indx)
stop
return
end
