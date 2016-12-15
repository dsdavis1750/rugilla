x = findgen(1000) + 1.
x = x / 10.
y1 = (x/2)*sqrt(x*x + 1.) - 0.5*(alog10(x + sqrt(x*x+1.)))
plot_oo,x,y1*5.43e-3
