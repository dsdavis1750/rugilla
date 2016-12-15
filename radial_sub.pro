pro radial_sub,subim,x_cen,y_cen,dim,rads,fav,x,line,coeffs,slope_err,minrange,maxrange,error,wing=wing

; Revised from Brian McNamara's prof3.pro.
; 01/21/99 by Alex Ware
;
;***********************************************************************
; this routine is for creating ellipse images and 
; annular profile maps of images
;***********************************************************************

;------------------------------------------------
lp=0  & def=-666.00 & badpix=def & npix=100
scl=1.   ;arcsec/pix
;____________________ INPUT DATA ______________________________________
gradfit=1 ; (N=0, Y=1) fit line to color profile                !!SET!!
prfplt=1  ; (N=0, Y=1) plot surface brightness profile          !!SET!!
bdpx=0    ; (N=0, Y=1) bad areas of chip removed                !!SET!!
xgal=[85.,0]  & ygal=[28.,0]; 			                !!SET!! 
eps=0.0 & phi=0.
th1=0. & th2=0.  ;opening slice between which to reject data     !!SET!!

wghts=[1,1,1,1,1,1, 1,1,1,1,1, 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]
ans=''

;------------------------------------------------------
szim=size(subim)
delimx=szim(1) & delimy=szim(2)	; x and y max of roi
;------------------------------------------------------
;******************************************************
; CM centroid - find image on chip
;******************************************************
;  margin vectors
xmv=total(subim,1)
ymv=total(subim,2)
szx=size(xmv) & szy=size(ymv)
xax=indgen(szx(1))
yax=indgen(szy(1))
mymv=total(xmv * xax)/total(xmv)
mxmv=total(ymv * yax)/total(ymv)

print,'radial_sub'
print,'CM Centroid Whole Image'
print,'Xcen=',mxmv,'Ycen=',mymv

xgal(0)=round(mxmv) & ygal(0)=round(mymv)
;----------------------------------
;if (bdpx eq 1) then badpixzw,subim
;-----------------------------------
; FIND SUBIMAGE CENTROID BY DAOPHOT CENTROIDER
cntrd,subim,mxmv,mymv,x_cen,y_cen,1
print,'DAOPHOT  CM centroid'
print,'Xcen=',x_cen,'Ycen=',y_cen

;******************************************************
; Change this depending on which centroid you want 
;******************************************************
y0=round(float(mxmv))  &  x0=round(float(mymv))
;y0=round(float(x_cen))  &  x0=round(float(y_cen))

;print,'Using Centroid, X=',y0,'Y=',x0

;******************************************************
; Cut out vertical stripe in wing images
;      (you may want to hardwire the xcen coords)
;******************************************************
;if keyword_set(wing) then begin
;	width = 4	; vertical stripe is 6 pix wide
;	xmin = xcen-(width/2) & xmax = xcen+(width/2)
;	ymin1 = 0 & ymax1 = 12
;	ymin2 = ymax1+13 & ymax2 = delimy-1
;	subim(xmin:xmax,ymin1:ymax1) = def
;	subim(xmin:xmax,ymin2:ymax2) = def
;endif
;------------------------------------------------------
; this xrcf image had a bad spot
;------------------------------------------------------
;if(trwid eq 'H-IAI-CR-1.003') then begin
;	xmin=(dim/2)-5 & xmax=(dim/2)+5 & ymin=1 & ymax=10
;	subim(xmin:xmax,ymin:ymax) = def
;endif
;------------------------------------------------------

;******************************************************
; Begin creating semimajor axis and theta maps
;******************************************************
profmsk=fltarr(delimx,delimy,2)
epspar=(1./(1.-eps))^2

  for i=0,delimx-1   do begin
        ydif=float(i-y0)
    for j=0,delimy-1 do begin
        xdif=float(j-x0) 
        
        if ((xdif eq 0) and (ydif eq 0)) then begin
           a=0 & thetd=0
           goto, elskp
         endif
       
         r=sqrt(xdif*xdif + ydif*ydif)
         numer=(ydif/r) & denom=(xdif/r)
         thet=atan(numer,denom)-phi*(!pi/180.) 
         asq=(epspar*sin(thet)^2 + cos(thet)^2)*r^2
         a=sqrt(asq)*scl  ;in arcsec
         thetd=thet*(180./!pi)
         if(thetd le 0.) then thetd= phi + thetd+180. ;360
         elskp:
         profmsk(i,j,0)=a & profmsk(i,j,1)=thetd
     endfor
   endfor

;******************************************************
; COMPUTE ANNULAR PROFILES
;******************************************************
rads=indgen(round(max(profmsk(*,*,0))))
;delrads=0.5
delrads=1.0
fsum=fltarr(n_elements(rads)) & fav=fltarr(n_elements(rads))
ftot=fltarr(n_elements(rads)) & nelements=fltarr(n_elements(rads))

  subim1=subim(*,*)

;******************************************************
; Make it start calculating radii at edge of central 
; stripe for wing images, otherwise it'll choke 
;******************************************************
if keyword_set(wing) then lowii=width/2 else lowii=0

for ii=lowii,n_elements(rads)-1 do begin

;PROFILOGIC &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& 
	favw=where((profmsk(*,*,0) ge rads(ii)-delrads)$
	and (profmsk(*,*,0) lt rads(ii)+delrads))
	if( (th1 ne 0.) and (th2 ne 0.))then begin 
		subout=where((profmsk(*,*,1) ge th1) and (profmsk(*,*,1) le th2))
		subim1(subout)=def
	endif

;&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&    
	subim2=subim1(favw)
	favw1=subim2(where(subim2 ne def))
	n0=n_elements(favw1)

	if(ii eq 0)then begin
		fav(ii)=favw1	; ALEX 6/4/99
	endif

	if(ii gt 0)then fav(ii)=avg(favw1)

	ftot(ii)=total(favw1)   
	nelements(ii)=n0
endfor

cumul=fltarr(n_elements(fav))
cumul(*)=0


for jjj=0, n_elements(fav)-1 do begin
if jjj eq 0 then  cumul(jjj)=total(ftot(jjj))
if jjj gt 0 then cumul(jjj)=cumul(jjj-1) + total(ftot(jjj))
 
endfor
cum_norm=cumul/cumul(n_elements(fav)-1)
;cum_norm=cumul/total(im)

rads=rads(where (ftot ne 0))
ftot=ftot(where ( ftot ne 0))
errorsp=(ftot + (ftot/sqrt(ftot)))/nelements
errorsm=(ftot - (ftot/sqrt(ftot)))/nelements
;rads=rads*24. ; in microns
;!fancy=3
;axis,yaxis=1,ystyle=1
;plot_oo,rads,fav,/ynozero,ystyle=8,yrange=[.0001,1000],xrange=[1,500],$
;  xtitle='!6R (pix)',ytitle='Cts/pix', title='HRMA-ACIS2C',/data

;SUBTRACT BACKGROUND

;print,'SUBTRACT BACKGROUND'
;!fancy=7

;fav=fav-0.232848
fav=fav-0.0

;******************************************************
; Fit a line to the outer section of the radial plot 
;  for wing data
; Adjust the minrange and maxrange params. to taste
;******************************************************
if keyword_set(wing) then begin
;	minrange=10 & maxrange=20
	minrange=10 & maxrange=15
;	minrange=15 & maxrange=20

	x=alog10(rads(minrange:maxrange))
	y=alog10(fav(minrange:maxrange))
	coeffs=poly_fit(x,y,1,line,errors,sigma) ;coeffs(0)=b , coeffs(1)=m

	res = y - line				;calc. residuals

	x = 10^x & y=10^y & line = 10^line	;convert to logspace
	errors = 10^errors & sigma = 10^sigma
	res = 10^res

	slope_err=stddev(res)
	print,'slope',coeffs(1)
	print,'slope error (logspace)',slope_err

;*************************************************************
; Plot the fit and residuals and all that - 
;  this may need to move to radial.pro
;*************************************************************

	!p.multi=[0,1,2]
	plot_oo,x,y,psym=8,xrange=[9,35],xstyle=1,$
	xtickname=[' ',' ',' ',' ',' ',' ',' '],ymargin=[0,3],$
	ytitle='avg. cts/pix',title='fit and residuals'
	oplot,x,line
	plot_oo,x,res,psym=4,ymargin=[3,0],xrange=[9,35],$
	xstyle=1,ytitle='avg. cts/pix',xtitle='radius (pix)'
	oplot,[1,40],[1,1],linestyle=1
	!p.multi=[0,1,1]

endif

totcts = total(ftot)

;*********************************************************
; Figure the errors
;*********************************************************
ymin=0.001
error=(sqrt(abs(ftot)))/nelements

;*********************************************************
skkk:

close,50
close,65
return
end

 
 
