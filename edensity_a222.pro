pro edensity,dummy
;
; physical constants
   c = 2.9979d5            ; km/s
   G = 6.67d-8
   k = 1.38d-16            ; erg/K
   Mpc = 3.08d24           ; cm
   Kpc = 3.08d21           ; cm
   M_sun = 1.989d33        ; grams
   proton_mass = 1.67d-24  ; grams
   erg_per_ev = 1.602d-12
   pi = 3.14159d0
   z = 0.2126
;
; Observational parameters
;
H     = 70   ; km/s/Mpc
scale = 3.26 ; arcsec/kpc
;
; fitted parameters
;
xspec_norms = [1.88251E-04,1.80440E-04,1.94876E-04,$
               2.07255E-04,1.83010E-04,2.16595E-04]
radii_inner = [0,40.,60.,80,100.,120.]
radii_outer = [40.,60.,80.,100.,120.,145.]
shell = dblarr(6)
;
; Compute the volumns of the shells
;

radii_inner = radii_inner*0.5 ; pixels to arcsec (Chandra)
radii_outer = radii_outer*0.5 ; pixels to arcsec (Chandra)

for i = 0,5 do begin
shell(i) = (4./3.)*pi*(radii_outer(i)^3 - radii_inner(i)^3)*scale^3*$
           Kpc^3
v = c * ((z+1)^2-1)/((z+1)^2+1)
D_a = v/H

e_dens2 = 1.d14*xspec_norms(i)*4*pi*(D_a*(1.+z))^2/shell(i)
n_dens = sqrt(e_dens2)/proton_mass
avg_rad = 0.5*(radii_inner(i)+radii_outer(i))

print,avg_rad*scale/1000.,n_dens
endfor
        
return
end
