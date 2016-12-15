pro n4261_mass,dummy
common pars,delta
x = findgen(1000) + 1.
x = x / 10.
gamma = 1.
beta = .311
delta = 3./2 * beta /( 1. - 0.2*(gamma - 1.))
print,delta
phi = delta * ( gamma - 1.)
a = .026
t10 = .85/10.
mu = 1.16
mean = 0.6
n0 = 6.77e-4
mvir = 1.92e14 * delta * gamma * t10 * (a/0.25)*(x^3)*(1+x^2)^(-(1.+phi))
mgas = mvir - mvir
mgas1 = 4.5e13*(a/0.4)^3*n0/.002 * (x/2)*sqrt(x*x + 1.) - 0.5*(alog10(x + sqrt(x*x+1.)))
!xtitle = 'Radius (Core radii)'
!ytitle = 'Mass (Solar Masses)'
!p.title = 'Mass profile for N4261 group'
!x.style=1
!x.range=[1.0,60]
;!y.range=[1.e8,2.e16]
!p.linestyle = 0
plot_oo,x,mvir
!p.linestyle = 1
;oplot,x,mgas1
for i = 0,999 do begin
int = simpson('density',0,x(i),count,TOL=.1)
mgas(i) = 3.11e17 * mu * n0 * a^3 * int
endfor
!p.linestyle = 2
oplot,x,mgas
!p.linestyle = 3
oplot,[23.81,23.81],[1.e5,1.e16]
print,mvir(0),mgas(0),mgas1(0)
print,mvir(238),mgas(238),mgas1(238)
;average density
avg_den = 9.64e-18*mgas(238)/(.622)^3
print,'Average number density',avg_den
!p.linestyle=1
oplot,[4.48,4.48],[1e11,1e14]
oplot,[10.24,10.24],[1e11,1e14]
oplot,[50.2 ,50.2],[1e11,1e14]
return
end


function density,y
common pars,delta
return, y*y/(1.+y*y)^delta
end
