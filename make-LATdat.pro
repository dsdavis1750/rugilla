;
;       MAKE-LATdat.PRO				Version:  OCTOBER 31, 2005
;
;       SubProgram				Calls, actions:
;
;	- setup					:sets up plotting device, colors, etc.
;	- randomBM				:generates repeatable, normally distributed random numbers
;	- photsrite				:writes detected photon log file
;	- Aeff					:generates effective area vs. energy, thin, thick, total
;	- Eresolution				;generates energy resolution vs. energy
;	- reconMode				;gets co-factors for angular and energy resolution
;	- disp_angles				;applies PSF per photon, thin or thick radiator
;	- disp_energy				;applies energy disperson per photon
;	- detect				;detects photons, applies Aeff(E,zenang), calls disp_angles, disp_energ
;	- makesymbol				;makes user defined symbols
;	- plotphots				;plots all the photons in cartesian coords
;	- sim_source				;does a source or bckgnd


pro setup
common mode, screen,dev,disply,dops,dosym,pause,ft,sz,axes,lulog
common colrtab,ltblue,blue,dkblue,green,blugreen,yellow,orange,red,black,cols
common lstyles, lnsty1,lnsty2,lnsty3,lnsty4,lnsty5,lnsty6,realthk,thicker,lnstys

;  set up colors, device, window, global plot parameters, font, linestyles

;getcolors, ncolors
black = 0  &  ltblue = 1  &  blue = 2  &  dkblue = 3  &  green = 4
blugreen = 5  &  yellow = 6  &  orange = 7  &  red = 8
cols = [blue,green,red,blugreen,black,orange,ltblue,dkblue,yellow]

if (dops eq 'n') then begin
;   set_plot, 'WIN'
   set_plot, 'X'
   xsizer = 850  &  ysizer = 850
   window, 0, title='Source+Bckgnd photons', retain=2,    $
      xpos=400,ypos=0,xsize=xsizer, ysize=ysizer, colors=ncolors
   !p.thick = 1  &  realthk = 2  &  thicker = 1
   sz = 1.0  &  !p.charthick = 1  &  !p.charsize = 1.2
   bckcol = orange
 endif else begin
   set_plot, 'ps'
   !p.position = [0.15,0.15,0.90,0.90]
   !p.thick = 2  &  realthk = 2  &  thicker = 2
   sz = 1.5  &  !p.charthick = 2  &  !p.charsize = 0.9
   bckcol = black
   black = 0  &  ltblue = 0  &  blue = 0  &  dkblue = 0  &  green = 0
   blugreen = 0  &  yellow = 0  &  orange = 0  &  red = 0
   cols = [blue,green,red,blugreen,black,orange,ltblue,dkblue,yellow]
 endelse

ft = '!6'                             ;Complex Roman

;  linestyles: 0: solid  1: ...  2: ---  3: -.-.  4:-...  5: ___ (long dash)

lnsty1 = 0  &  lnsty2 = 0  &  lnsty3 = 0  &  lnsty4 = 0
lnsty5 = 0  &  lnsty6 = 0  &  lnsty7 = 0

return  &  end



function randomBM, iseed, Nget

;  see http://en.wikipedia.org/wiki/Box-Muller_transform

Ngen = ceil(Nget/2.)

R   = randomu(iseed, Ngen)
phi = randomu(iseed, Ngen)

; (unlike IDL's randomn) z0 and z1 are repeatable, given the same seed

z0 = cos(2.*!PI*phi) * sqrt(-2. * alog(R))
z1 = sin(2.*!PI*phi) * sqrt(-2. * alog(R))

zarr = [z0, z1]
zarr = zarr(0:Nget-1)
if (N_elements(zarr) eq 1) then zarr = zarr(0)

return, zarr
end



pro photsrite
common mode, screen,dev,disply,dops,dosym,pause,ft,sz,axes,lulog
common angles, zenang, aziang
common photons, Energy, Angle, phottype
common sources, Npntsour, Nsours, betas, zenangs, aziangs

;  write log file for detected LAT burst

Ndet = n_elements(phottype)

;prefix = 'D:\idlprog\Glast\SourceFinder\data\'
prefix = './data/'
logfilesp = prefix + 'LAT-Nsour-bck.lis'

openw, lulog, logfilesp, /get_lun

printf, lulog, f="('  Number of point sources = ',I3)", Npntsour
printf, lulog, f="('   Betas   AziAngs  Zenangs ')"
for i = 0,Npntsour-1 do printf, lulog, f="(3F9.4)", betas(i), zenangs(i), aziangs(i)

printf, lulog, f="('  Nphots Detect  = ',I9)", Ndet
printf, lulog, f="('      Energy      Zenith       Azimuth   Type:')"

matrix = fltarr(4, Ndet)
matrix(0,*) = Energy
matrix(1,*) = Angle(*,0)
matrix(2,*) = Angle(*,1)
matrix(3,*) = phottype
printf, lulog, f="(3F13.6,3X,I2)", matrix

close, lulog  &  free_lun, lulog




if (pause eq 'y') then begin
   ques = ''  &  read, '? ', ques
   if (strlowcase(strmid(ques,0,1)) eq 's') then stop, 'examine matrix info, in photsrite'
 endif

return  &  end



pro Aeff, nptsE
common LATcharacter, logE, Aefftot, Aeffthin, Aeffthick, Eres
common ModeOption, Onboard, Gleam, factPSF, factEres
common seeder, iseed

; make effective area arrays, in cm^2, as function of energy; then convert
; to detection probabilities given illumination area of 6 m^2 in GRBfullsim_15.

; Rectangle Grid:
;  X,   Y : logE, cm^2
;197, 686 = 1.0,     0
;197, 108 = 1.0, 20000
;899, 686 = 5.0,     0
;899, 108 = 5.0, 20000

X0 = 197.
Y0 = 686.
deltaX = 899. - 197.
deltaY = 108. - 686.

logE0 = 1.0
area0 = 0.0
deltlogE = 5.0 - 1.0
deltarea = 2E4 - 0.0

logEfact = deltlogE / deltaX
afact    = deltarea / deltaY

decs = 4.  &  nperdec = 25
nptsE = float(round(decs * nperdec)) + 1.
logE = fltarr(nptsE)

Aeffthin  = fltarr(nptsE)
Aeffthick = fltarr(nptsE)
Aefftot   = fltarr(nptsE)

; general formulae:
; logE = logE0 + Efact*deltaXlin
; area = area0 + afact*deltaYlin

; piece-wise continuous function for thin radiator section
; line segments:
; (X1,Y1) -> (X2,Y2)

line = fltarr(16,4)
line( 0,*) = [197, 686,  210, 680]
line( 1,*) = [210, 680,  253, 659]
line( 2,*) = [253, 659,  310, 630]
line( 3,*) = [310, 630,  392, 591]
line( 4,*) = [392, 591,  433, 574]
line( 5,*) = [433, 574,  476, 558]
line( 6,*) = [476, 558,  519, 542]
line( 7,*) = [519, 542,  564, 529]
line( 8,*) = [564, 529,  609, 520]
line( 9,*) = [609, 520,  655, 512]
line(10,*) = [655, 512,  700, 505]
line(11,*) = [700, 505,  745, 501]
line(12,*) = [745, 501,  792, 500]
line(13,*) = [792, 500,  837, 497]
line(14,*) = [837, 497,  899, 494]
line(15,*) = [899, 494,  951, 494]

Xrange = line(14,2) - line(0,0)

for ipnt = 0,nptsE-1 do begin
   dx = (ipnt/(nptsE-1)) * Xrange
   logE(ipnt) = logE0 + logEfact * dx

   Xpnt = X0 + dx
   dum = where( (Xpnt ge line(*,0)) AND (Xpnt lt line(*,2)), ndum)
   m = dum(0)
   slope = (line(m,3) - line(m,1)) / (line(m,2) - line(m,0))
   Xoffset = Xpnt - line(m,0)
   Yoffset = slope * Xoffset
   Ypnt = line(m,1) + Yoffset

   Aeffthin(ipnt) = area0 + afact * (Ypnt - Y0)
 endfor

; piece-wise continuous function for thin+thick=total radiator sections

line = fltarr(16,4)
line( 0,*) = [197, 686,  210, 680]
line( 1,*) = [210, 680,  253, 646]
line( 2,*) = [253, 646,  310, 573]
line( 3,*) = [310, 573,  392, 476]
line( 4,*) = [392, 476,  433, 440]
line( 5,*) = [433, 440,  476, 411]
line( 6,*) = [476, 411,  519, 388]
line( 7,*) = [519, 388,  564, 368]
line( 8,*) = [564, 368,  609, 352]
line( 9,*) = [609, 352,  655, 338]
line(10,*) = [655, 338,  700, 326]
line(11,*) = [700, 236,  745, 318]
line(12,*) = [745, 318,  792, 311]
line(13,*) = [792, 311,  837, 308]
line(14,*) = [837, 308,  899, 306]
line(15,*) = [899, 306,  951, 306]

for ipnt = 0,nptsE-1 do begin
   dx = (ipnt/(nptsE-1)) * Xrange

   Xpnt = X0 + dx
   dum = where( (Xpnt ge line(*,0)) AND (Xpnt lt line(*,2)), ndum)
   m = dum(0)
   slope = (line(m,3) - line(m,1)) / (line(m,2) - line(m,0))
   Xoffset = Xpnt - line(m,0)
   Yoffset = slope * Xoffset
   Ypnt = line(m,1) + Yoffset

   Aefftot(ipnt) = area0 + afact * (Ypnt - Y0)
 endfor

Aeffthick = Aefftot - Aeffthin

; following lines used for check:
;Aeffthick = fltarr(101) + 6500.
;Aeffthin = Aeffthick
;Aefftot = Aeffthin + Aeffthick
;stop, 'examine Aefftot'

; If input # of photons was made for Gleam, then:
; convert area(E) into probability for detection:
; 6 m^2 disk was illuminated in GRBfullsim_10;
; else, normalize by frontal area = 160^2 cm^2

if (Gleam eq 'y') then begin
   Adisk = 60000.
   Aeffthin  = Aeffthin  / Adisk
   Aeffthick = Aeffthick / Adisk
   Aefftot   = Aefftot   / Adisk
 endif else begin
   Afront = 160. * 160.   ; = 25600.
   Aeffthin  = Aeffthin  / Afront
   Aeffthick = Aeffthick / Afront
   Aefftot   = Aefftot   / Afront
 endelse

;stop, 'examine Aefftot, Aeffthin, Aeffthick'

return  &  end



pro Eresolution, nptsE
common LATcharacter, logE, Aefftot, Aeffthin, Aeffthick, Eres
common seeder, iseed

; make energy resolution array, sigmaE/E, as function of energy

; Rectangle Grid:
;  X,   Y : logE,sigE/E
;240, 597 = 1.0,     0
;240,  47 = 1.0,   0.5
;891, 597 = 5.0,     0
;891,  47 = 5.0,   0.5

X0 = 240.
Y0 = 597.
deltaX = 891. - 240.
deltaY =  47. - 597.

logE0 = 1.0
sigE0 = 0.0
deltlogE = 5.0 - 1.0
deltsigE = 0.5 - 0.0

logEfact = deltlogE / deltaX
sfact    = deltsigE / deltaY

Eres  = fltarr(nptsE)

; piece-wise continuous function for "standard deviation" of energy resolution
; line segments:
; (X1,Y1) -> (X2,Y2)

line = fltarr(12,4)

line( 0,*) = [240, 317,  275, 372]
line( 1,*) = [275, 372,  333, 454]
line( 2,*) = [333, 454,  344, 466]
line( 3,*) = [344, 466,  368, 484]
line( 4,*) = [368, 484,  417, 504]
line( 5,*) = [417, 504,  482, 518]
line( 6,*) = [482, 518,  534, 525]
line( 7,*) = [534, 525,  598, 529]
line( 8,*) = [598, 529,  673, 522]
line( 9,*) = [673, 522,  780, 505]
line(10,*) = [780, 505,  891, 483]
line(11,*) = [891, 483,  944, 472]

Xrange = line(10,2) - line(0,0)

for ipnt = 0,nptsE-1 do begin
   dx = (ipnt/(nptsE-1)) * Xrange

   Xpnt = X0 + dx
   dum = where( (Xpnt ge line(*,0)) AND (Xpnt lt line(*,2)), ndum)
   m = dum(0)
   slope = (line(m,3) - line(m,1)) / (line(m,2) - line(m,0))
   Xoffset = Xpnt - line(m,0)
   Yoffset = slope * Xoffset
   Ypnt = line(m,1) + Yoffset

   Eres(ipnt) = sigE0 + sfact * (Ypnt - Y0)
 endfor

return

; *********************************

nphots = 1000000L
E = fltarr(nphots) + 2.

sigma2 = 2.*0.06*E(0)
sigL = 0.625*sigma2  &  sigH = 0.375*sigma2

lofrac = (sigL/(sigL+sigH))
LH = randomu(iseed, nphots)
dumL = where (LH le lofrac, ndumL)
dumH = where (LH gt lofrac, ndumH)

ErespL = E(dumL) - sigL*abs(randomBM(iseed,ndumL))
ErespH = E(dumH) + sigH*abs(randomBM(iseed,ndumH))
Eresp = [ErespL, ErespH]
Eresp = Eresp(sort(Eresp))

biner = 0.01  &  Emin = 0.  &  Emax = 3.
nbins = (Emax - Emin) / biner + 1
xhist = histogram(Eresp, min=Emin, max=Emax, binsiz=biner)
Eplot = biner*findgen(nbins)
!x.style = 1
plot, Eplot+biner/2, xhist, psym=10, xrange=[0.8,2.4]

stop, 'examine stuff'

return  &  end



pro disp_angles, indx, iphot, AngleDet
common angles, zenang, aziang
common LAT, logEinc, LATenergies
common LATcharacter, logE, Aefftot, Aeffthin, Aeffthick, Eres
common ModeOption, Onboard, Gleam, factPSF, factEres
common seeder, iseed

; disperse the incident angle, according to naive, canned PSFs for thick and thin radiator zones

dip = 180./!PI

Athinrat = (Aeffthin(indx)>1.E-8) / (Aefftot(indx)>1.E-8)

chanceradzone = randomu(iseed)
if (chanceradzone lt Athinrat) then begin              ;determine which layer type, thick or thin
   Err1sig =      factPSF * (10^(-0.95*(logEinc(iphot) + 1.523) + 1.0) + 0.03)
 endif else begin
   Err1sig = 1.35*factPSF * (10^(-0.95*(logEinc(iphot) + 1.523) + 1.0) + 0.03)
 endelse

; compute direction for detected photon, based on the single-Gaussian characterization of the PSF

rho = abs(Err1sig*randomBM(iseed,1)) / dip             ;distance of det'd photon from incident position
alpha = 2*!PI*randomu(iseed)                           ;azimuth of det'd position about incident position, rel. instr. axis
zinc = zenang / dip

CosZdet = cos(rho)*cos(zinc) + sin(rho)*sin(zinc)*cos(alpha)
zdet = Acos(CosZdet)
AngleDet(iphot,0) = dip * zdet                         ;detected zenith angle

beta = dip * Asin(sin(rho)*sin(alpha)/sin(zdet))       ;angle between zinc and zdet
AngleDet(iphot,1) = beta + aziang

;zsinfact = sin(zenang/dip)
;!x.style = 1  &  !y.style = 1
;plot, AngleDet(ndet(isim),1)+fltarr(1), zdet*dip+fltarr(1), psym=1, symsiz=0.125,    $
;      xrange=[aziang-5./zsinfact,aziang+5./zsinfact], yrange=[zenang-5.,zenang+5]
;!noeras = 1

return  &  end



pro disp_energy, indx, iphot, EnergyDet
common LATcharacter, logE, Aefftot, Aeffthin, Aeffthick, Eres
common LAT, logEinc, LATenergies
common seeder, iseed
common ModeOption, Onboard, Gleam, factPSF, factEres

; disperse the incident energy (Einc), using two-sided Gaussian, and sigma = Eres(Einc)

Einc = LATenergies(iphot)

sigma2 = 2.*Eres(indx)*Einc
sigL = 0.625*sigma2
sigH = 0.375*sigma2

lofrac = (sigL/(sigL+sigH))
LH = randomu(iseed)

if (LH le lofrac) then Eresp = Einc - factEres * sigL*abs(randomBM(iseed,1))   $
                  else Eresp = Einc + factEres * sigH*abs(randomBM(iseed,1))

EnergyDet(iphot) = Eresp > 0.010

return  &  end



pro detect, nptsE, Nphots, AngleDet, EnergyDet, hit
common angles, zenang, aziang
common LATcharacter, logE, Aefftot, Aeffthin, Aeffthick, Eres
common LAT, logEinc, LATenergies
common seeder, iseed

; per burst:  modify the Aeff, according to angle of incidence formula generated by Toby

dip = 180./!PI
Zenang1 = zenang / dip
cosZenang = cos(Zenang1)
a = 0.32  &  b = 0.80  &  c = 0.18

if (cosZenang gt b)    then FactAeff = (cosZenang + a) / (1 + a)                        else  $
   if (cosZenang gt c) then FactAeff = (cosZenang - c) * (b + a) / ((1 + a) * (b - c))  else  $
                            FactAeff = 0.

for iphot = 0,Nphots-1 do begin
   dum = where(logEinc(iphot) lt logE, ndum)
   if (ndum gt 0) then begin
      dum_1 = (dum(0) - 1) > 0
      chanceindx = randomu(iseed)
      if (chanceindx gt 0.5) then indx = dum(0)  else indx = dum_1
    endif else begin
      indx = nptsE-1
    endelse

   chancephot = randomu(iseed)
   if (chancephot lt FactAeff*Aefftot(indx)) then begin  ;photon detected, so apply dispersions
      hit(iphot) = 1                                     ;this photon detected
      disp_angles, indx, iphot, AngleDet                 ;apply PSFs, thick and thin radiator zones
      disp_energy, indx, iphot, EnergyDet                ;apply energy dispersion to incident energy
    endif

 endfor

return  &  end



pro make_LATenergies, beta, Nphots
common LAT, logEinc, LATenergies
common seeder, iseed

; make the energies according to a power-law spectrum

Ethres = 0.02   ;minomum energy detected
Emax   = 100.   ;maximum energy allowed = 100 GeV

if (beta eq 1.0) then begin         ;handle the one special case
   LATenergies = Ethres * exp(randomu(iseed,Nphots) * alog(Emax/Ethres))
 endif else begin
   X = 1. - exp((1. - beta) * alog(Emax/Ethres))
   LATenergies = Ethres * exp(alog(1. - X*randomu(iseed,Nphots))/(1. - beta))
 endelse

return  &  end



pro makesymbol, iplt, symsiz, vertices

;  make one of five polygons:  diamond, triangle, square, circle, upsidedown triangle

MODer = iplt MOD 5

symsizes = fltarr(5) + 1.5

Nvertices = [5, 4, 5, 33, 4]
Nphases = [0., 90., 45., 0., 270.]

vertices = findgen(Nvertices(MODer)) * (!PI * 2. / (Nvertices(MODer)-1))   $
         + Nphases(MODer) * !PI/180.

symsiz = symsizes(MODer)

return  &  end



pro plotphots, dip, scale
common mode, screen,dev,disply,dops,dosym,pause,ft,sz,axes,lulog
common colrtab,ltblue,blue,dkblue,green,blugreen,yellow,orange,red,black,cols
common lstyles, lnsty1,lnsty2,lnsty3,lnsty4,lnsty5,lnsty6,realthk,thicker,lnstys
common photons, Energy, Angle, phottype

makesymbol, 3, symsiz, vertices
usersym, cos(vertices), sin(vertices), fill=1

icols = [yellow, orange, red]

!p.title = ft
!x.title = ft + 'Xcart'
!y.title = ft + 'Ycart'

if (scale eq 's') then begin
   !y.range = [-0.15,0.35]  &  !y.ticks = 10
   !x.range = [-0.25,0.25]  &  !x.ticks = 10
 endif else begin
   !y.range = [-1,1]  &  !y.ticks = 10
   !x.range = [-1,1]  &  !x.ticks = 10
 endelse

!x.style = 1  &  !y.style = 1
!p.position = [0.13,0.13,0.95,0.95]

dumpltS = where(phottype gt 0, Ns)
dumpltB = where(phottype lt 0, Nb)

Xcart = sin(Angle(*,0)/dip) * cos(Angle(*,1)/dip)
Ycart = sin(Angle(*,0)/dip) * sin(Angle(*,1)/dip)
logE = alog10(Energy)
Esiz = 0.75/(2.5 + logE) > 0.


npts = n_elements(Xcart)
openw, 66, 'LAT-data_cart.txt'
for i = 0,npts-1 do printf, 66, f="(2F9.4)",  Xcart(i)+1., Ycart(i)+1.
close, 66  &  free_lun, 66



plot, /nodata, fltarr(1), fltarr(1), col=orange
for i = 0,Ns-1 do oplot, Xcart(dumpltS(i))+fltarr(1), Ycart(dumpltS(i))+fltarr(1), psym=8, symsiz=Esiz(i),    $
                  col=icols(phottype(dumpltS(i))-1)
for i = 0,Nb-1 do oplot, Xcart(dumpltB(i))+fltarr(1), Ycart(dumpltB(i))+fltarr(1), psym=8, symsiz=Esiz(i), col=blue

if (pause eq 'y') then begin
   ques = ''  &  read, '? ', ques
   if (strlowcase(strmid(ques,0,1)) eq 's') then stop, 'examine plot'
 endif

return  &  end



pro sim_source, Nsim, beta, EnergyDet1, AngleDet1, typephot1, BorS, nptsE, Nsimdet
common LAT, logEinc, LATenergies
common angles, zenang, aziang

AngleDet =  fltarr(Nsim,2)
EnergyDet = fltarr(Nsim)
hit       = intarr(Nsim)                                  ;detect (1) or not (0)

make_LATenergies, beta, Nsim                              ;makes incident source photons
logEinc = alog10(LATenergies)                             ;array of log(incident energies)
detect, nptsE, Nsim, AngleDet, EnergyDet, hit             ;detect photons

dumhits = where(hit eq 1, ndum)
Nsimdet = ndum

if (ndum gt 0) then begin
   if (BorS eq 'S') then begin
      print, 'Number of source photons detected = ', ndum
      typephot1  = fltarr(Nsimdet) + 1
    endif else begin
      print, 'Number of bckgnd photons detected = ', ndum
      typephot1  = fltarr(Nsimdet) - 1
    endelse

   EnergyDet1 = EnergyDet(dumhits)
   AngleDet1  = AngleDet(dumhits,*)
 endif else begin
   print, 'Zero source hits:  abort'
   stop, ''
 endelse

return  &  end


;  start Main Program
;
print, f="(/,' Simulates detection of GLAST LAT sources, based on characterizations of LAT response',/)"
common mode, screen,dev,disply,dops,dosym,pause,ft,sz,axes,lulog
common colrtab,ltblue,blue,dkblue,green,blugreen,yellow,orange,red,black,cols
common lstyles, lnsty1,lnsty2,lnsty3,lnsty4,lnsty5,lnsty6,realthk,thicker,lnstys
common plotpars,ptit,xtit,ytit,pltpnts,symtype,xmaxer,ysym
common LAT, logEinc, LATenergies
common angles, zenang, aziang
common LATcharacter, logE, Aefftot, Aeffthin, Aeffthick, Eres
common seeder, iseed
common ModeOption, Onboard, Gleam, factPSF, factEres
common photons, Energy, Angle, phottype
common sources, Npntsour, Nsours, betas, zenangs, aziangs

;iseed = 92962880521L
iseed = 962880521L
dip = 180./!PI

dev = 'WIN'  &  disply = 'c'  &  screen = 'y'  &  dops = 'n'

parseme = 'nysy'  &  ques = ''
print, '                      Default parseme = ', parseme
read, '              Psplt Pause Scale WrtLog: ', ques
if (ques ne '') then parseme = strlowcase(ques)
dops = strmid(parseme,0,1)   &  pause = strmid(parseme,1,1)
Scale = strmid(parseme,2,1)  &  wrtlog = strmid(parseme,3,1)

setup

Gleam = 'n'

Aeff, nptsE                                                ;get run of effective areas vs. energy
Eresolution, nptsE                                         ;get run of energy resolution vs. energy

factPSF = 1.
factEres = 1.
logE = logE - 3.                                           ;convert: MeV to GeV

;************************************************************************************************************
; simulate bckgnd photons, whose spatial distribution will be flat;
;   but whose energy distribution reflects LAT's Aeff(E) and bckgnd spectrum.
; so, same sequence and machinery as for source photons, except the bckgnd
;   photons' positions will be randomized after the fact.
;************************************************************************************************************

Nbck = 6000
beta  = 2.0

zenang = 10.
aziang = 75.

sim_source, Nbck, beta, EnergyDet_B, AngleDet_B, typephot_B, 'B', nptsE, Nbckdet

; randomize the positions of the bckgnd photons

azirand = 360.*randomu(iseed, Nbckdet)
zenmax = 50.
Zen_Norm = cos(0.) - cos(zenmax/dip)
cos_Zen = 1. - randomu(iseed, Nbckdet) * Zen_Norm
zenrand = acos(cos_Zen) * dip
AngleDet_B(*,0) = zenrand
AngleDet_B(*,1) = azirand

; put source and bckgnd photons into one array, then sort them by energy

Energy   = EnergyDet_B
Angle    = AngleDet_B
phottype = typephot_B

;************************************************************************************************************
; simulate detected point source photons, whose spatial and energy dists
;   reflect LAT's Aeff(E) and PSF(E), as well as the source's spectrum.
;************************************************************************************************************

Npntsour = 3

Nsours  = [600, 600, 600]
betas   = [2.0, 2.0, 2.0]
zenangs = [10., 10., 10.]
aziangs = [75.,  0.,150.]


for i = 0,Npntsour-1 do begin
   zenang = zenangs(i)  &  aziang = aziangs(i)
   sim_source, Nsours(i), betas(i), EnergyDet_S, AngleDet_S, typephot_S, 'S', nptsE, Nsourdet

   Energy   = [Energy,   EnergyDet_S]
   Angle    = [Angle,    AngleDet_S]
   phottype = [phottype, typephot_S+i]
 endfor

sorter   = sort(Energy)
Energy   = Energy(sorter)
Angle    = Angle(sorter,*)
phottype = phottype(sorter)

;************************************************************************************************************
; plot time:  plot the source and bckgnd photons in a cartesian coord system
;************************************************************************************************************

plotphots, dip, scale

; write out the file of detected photons:  angles and energies

if (wrtlog eq 'y') then photsrite

thender:
end
