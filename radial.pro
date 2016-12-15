pro radial,rois,xcenters,ycenters,dim,files,rads_arr,fav_arr,err_arr,wing=wing

; Revised from Brian McNamara's prof3.pro.
; 3/99 by Alex Ware
;
; this routine is for creating ellipse images and 
; annular profile maps of images
;

xdim=dim & ydim=dim

rads_arr=fltarr((n_elements(files)),120)
fav_arr=fltarr((n_elements(files)),120)
err_arr=fltarr((n_elements(files)),120)

for j=0,(n_elements(files)-1) do begin

roi=rois[*,*,j]
filename=files[j]
xcen=xcenters[j]
ycen=ycenters[j]

;------------------------------------------------------
; call the radial_sub program
;------------------------------------------------------
if keyword_set(wing) then begin
radial_sub,roi,xcen,ycen,dim,rads,fav,x,line,coeffs,slope_err,minrange,maxrange,error,/wing
endif else begin
radial_sub,roi,xcen,ycen,dim,rads,fav,x,line,coeffs,slope_err,minrange,maxrange,error
endelse

rads_arr[j,0:(n_elements(rads)-1)]=rads
fav_arr[j,0:(n_elements(fav)-1)]=fav
err_arr[j,0:(n_elements(error)-1)]=error

;if keyword_set(wing) then begin
;need to output: x,line,coeffs,slope_err,minrange,maxrange


skkk:

endfor

end

 
 







