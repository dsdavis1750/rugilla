pro acisi_psf_interpol, energy,roff,r,psf
; read the psfsize table currently /home/ciao/data/psfsize20010416.fits
; and linearly interpolate the psf vs r from the table. 
;
; Yuxuan Yang 04/21/03
;
;PSF FILE DIR is set to the default ciao dir
;
;psftable='/astromake/opt/ciao/3.0.2/data/psfsize20010416.fits'
psftable='/usr1/local/Chandra/ciao_3.0/data/psfsize20010416.fits'
;
;In the PSF table there are 9 off-axis angles and 5 energy bands
theta=[0,1,2,3,5,7,10,15,20]
E=[0.277,1.4967,4.51084,6.40384,8.63886]
;
;
r=fltarr(5,9,495)
frac=fltarr(5,9,495)
FITS_OPEN,psftable,fcb
;
; First index is energy
; Second index is off-axis angle
;
for i = 0,4 do begin
	for j = 0,8 do begin
	extno=i*9+j+6
	FITS_READ,fcb,tab,h,exten_no=extno
	radius=TBGET(h,tab,'RADIUS')
	fraction=TBGET(h,tab,'FRACTION')
	r(i,j,*)=radius
	frac(i,j,*)=fraction
;	print,extno,E(i),theta(j)
	endfor
endfor	
FITS_CLOSE,fcb
;Interpolate the PSF size file to the desired energy and 	
;1) interpol to energy
	r1=fltarr(9,495)
	frac1=fltarr(9,495)
	for i = 0, 8 do begin 			;For each theta
		for j = 0, 494 do begin		;For each radius
			result=interpol(frac(*,i,j),E,energy)
			r1(i,j)=r(0,i,j)
			frac1(i,j)=result
		endfor
	endfor
;
;2) Interpol to roff
;
	r2=fltarr(495)
	frac2=fltarr(495)
	for i = 0, 494 do begin			;For each radius
		result = interpol(frac1(*,i),theta,roff)
		r2(i) = r1(0,i)
		frac2(i) = result
	endfor
r = r2
psf = frac2
return
end
		  
