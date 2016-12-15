pro mass,dummy
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
rcore = 1.303e-3  ; Mpc
rmax = 3.0e-2     ; Mpc
temp = 0.746      ; keV
mu = 0.62         ; mean molecular weight for fully ionized gas
beta = 0.486      ; slope of surface brightness
n0 = .332         ; central electron number density

;
; radial vector out to 100
x = findgen(1000) + 1.
x = x / 1000.       ; set range from 0 - 1. Mpc
x = x /rcore        ; convert x to units of rcore
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
r = x * Mpc
kT = (temp * 1000.0) * erg_per_ev
sigma_sq = kT/(mu*proton_mass)
x_factor = (x^3) / (1 + x^2)
m_tot = (3*beta*sigma_sq*x_factor*rcore*Mpc)/G - (x^2 *k* sigma_sq/(G*kT))*dtdx 
;m_tot = (-1.*k*x*x*Mpc/(G*mu*proton_mass))*(dtdr*rcore^2 - 3*beta*temp*(x*rcore/(1. + x*x)))
;m_tot = 1.3e7*m_tot / M_sun ; convert to solar masses
m_tot = m_tot / M_sun ; convert to solar masses
;
; electron distribution
n_elec = n0 * ( 1. / (1. + x^2))^(-1.5 * beta)
;
entropy = alog10(kT/k / n_elec^(2./3.))
;
!p.multi=[0,0,2]
!xtitle = 'Radius (kpc)'
!ytitle = 'Mass (Solar Masses)'
!p.title = 'Mass profile for NGC 1407'
!p.linestyle = 0
plot_oo,x*rcore*1000.,m_tot,yrange=[1.e7,1.e14],xrange=[1.0,1000.]
!p.linestyle = 8
oplot,x*rcore*1000.,mgas
!p.linestyle = 0
plot_oi,x*rcore*1000.,entropy
!p.multi=0
print,entropy
return
end
