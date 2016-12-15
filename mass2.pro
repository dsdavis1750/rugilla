pro mass2,dummy
x = findgen(1000) + 1.
x = x / 10.
gamma = 1.
beta = .84
delta = 3./2 * beta /( 1. - 0.2*(gamma - 1.))
phi = delta * ( gamma - 1.)
a = .595
t10 = .8
n0 = 1.75e-3
!xtitle = 'Radius (Core radii)'
!ytitle = 'Mass (Solar Masses)'
!p.title = 'Mass profile for A400'
mvir = 1.92e14 * delta * gamma * t10 * (a/0.25)*(x^3)*(1+x^2)^(-(1.+phi))
mgas = 4.5e13*(a/0.4)^3*n0/.002 * ( x - atan(x))
!p.linestyle = 0
plot_oo,x,mvir
!p.linestyle = 8
oplot,x,mgas
return
end
