pro eer,rois,xcenters,ycenters,dim,files,radial_arr,ee_arr,eer_arr,wing=wing

; Outputs the encircled energy for an input image array
; Last revised: 5/99 by Alex Ware
; Property of the Smithsonian Astrophysical Observatory

xdim=dim & ydim=dim

radial_arr=fltarr((n_elements(files)),100)
ee_arr=fltarr((n_elements(files)),100)
eer_arr=fltarr((n_elements(files)),4)

for j=0,(n_elements(files)-1) do begin

roi=rois[*,*,j]
filename=files[j]
xcen=xcenters[j]
ycen=ycenters[j]

;------------------------------------------------------
find_max, roi, xcen, ycen

szim=size(roi)
delimx=szim(1) & delimy=szim(2)	; x and y max of roi
;------------------------------------------------------
; The PSF wing images tend to have a vertical stripe;
;  this section blocks it out
;------------------------------------------------------
;def=-666.00
def=0
if keyword_set(wing) then begin
	width = 6	; vertical stripe is 4 pix wide
;	xmin = xcen-(width/2) & xmax = xcen+(width/2)
	xmin = (xdim/2)-(width/2) & xmax = (xdim/2)+(width/2)
	ymin = 0 & ymax = delimy-1
	roi(xmin:xmax,ymin:ymax) = def
	subim=roi
endif
;------------------------------------------------------------------------
;Call the program that actually calculates the encircled energies

calc_encircled_energy, roi, xcen, ycen, radial_distance, c_sum, ee, error,nr

eer_arr[j,3]=nr		;Normalization radius

;------------------------------------------------------------------------
; put data into arrays for output to plotting subroutine
;-------------------------------------------------------
radial_arr[j,0:(n_elements(radial_distance)-1)]=radial_distance
ee_arr[j,0:(n_elements(ee)-1)]=ee

;**************************************************************
;Find 50%, 70%, 90% energy radii, for scaled and unscaled EE
;**************************************************************
l50=max(where(ee lt 0.5)) & h50=min(where(ee gt 0.5))
l70=max(where(ee lt 0.7)) & h70=min(where(ee gt 0.7))
l90=max(where(ee lt 0.9)) & h90=min(where(ee gt 0.9))

lo50 = ee(l50)
hi50 = ee(h50)
m50 = hi50-lo50		;slope
b50 = lo50-(m50*l50)	;intercept
e50 = (0.5-b50)/m50	;xvalue at y=0.5

lo70 = ee(l70)
hi70 = ee(h70)
m70 = hi70-lo70		;slope
b70 = lo70-(m70*l70)	;intercept
e70 = (0.7-b70)/m70	;x at y=0.7

lo90 = ee(l90)
hi90 = ee(h90)
m90 = hi90-lo90		;slope
b90 = lo90-(m90*l90)	;intercept
e90 = (0.9-b90)/m90	;x at y=0.9

;**************************************************************
;Calculate errors for encircled energy radii
;**************************************************************
err = error*max(c_sum)		; c_sum is cumulative sum (total cts)

y50=0.5*max(c_sum) & y70=0.7*max(c_sum) & y90=0.9*max(c_sum)
err50=sqrt(y50) & err70=sqrt(y70) & err90=sqrt(y90)

h50 = y50 + err50
m = m50*max(c_sum) & b = b50*max(c_sum)
e50_h = (h50-b)/m
e50_err = e50_h - e50

h70 = y70 + err70
m = m70*max(c_sum) & b = b70*max(c_sum)
e70_h = (h70-b)/m
e70_err = e70_h - e70

h90 = y90 + err90
m = m90*max(c_sum) & b = b90*max(c_sum)
e90_h = (h90-b)/m
e90_err = e90_h - e90

eer_arr[j,0]=e50
eer_arr[j,1]=e70
eer_arr[j,2]=e90

;-----------------------------------------------------
;Write encircled energy to output files 
;-----------------------------------------------------
;all grades: filename_all.dat
;std grades: filename_std.dat

printf,1,filename
printf,1,'e50',e50,' +-',e50_err
printf,1,'e70',e70,' +-',e70_err
printf,1,'e90',e90,' +-',e90_err

;**************************************************************


endfor
return
end
