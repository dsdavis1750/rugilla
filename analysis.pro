pro analysis,file_array,spec_array,spec_range,type_array,multi,outfile,std=std

; read in file, select region of interest,
; make a set of plots of analysis products,
; dump output numbers to text files
;
; written by Alex Ware 3/99
; Revised 6/99

;****************************************************************
; Set some things up
;****************************************************************
dim=100  ; dimensions of ROI

pw=strpos(file_array[0],'PW')
if(pw ne -1)then wing=1 else wing=0

roi_array=fltarr(dim+1,dim+1,n_elements(file_array))
xcen_array=fltarr(n_elements(file_array))
ycen_array=xcen_array

;****************************************************************
; Read in images, make ROIs from images. 
;  Can handle XRCF, flight, marx, or level 1 data.
;****************************************************************

for j=0,(n_elements(file_array)-1) do begin

type=type_array[j]
filename=file_array[j]

if(type eq 'xrcf')or(type eq 'flight')then begin
	file_all=filename+'_all.fits_img.fits'
	image=readfits(file_all)
	mkroi,image,dim+1,xcen,ycen,roi
	roi_array[0:(n_elements(roi[*,0])-1),0:(n_elements(roi[0,*])-1),j]=roi
	xcen_array[j]=xcen & ycen_array[j]=ycen
	image=0

;	Standard grades
	if keyword_set(std) then begin
	file_std=filename+'_02346.fits_img.fits'
	image=readfits(file_std)
	mkroi,image,dim+1,xcen,ycen,roi
	roi_array[0:(n_elements(roi[*,0])-1),0:(n_elements(roi[0,*])-1),j]=roi
	xcen_array[j]=xcen & ycen_array[j]=ycen
	image=0
	endif

endif

if(type eq 'marx')then begin

	grade=read_marx_file(filename+'grade.dat')
	pha=read_marx_file(filename+'pha.dat')
	x=read_marx_file(filename+'chipx.dat')
	y=read_marx_file(filename+'chipy.dat')
	image=make_image(x,y)
	mkroi,image,dim+1,xcen,ycen,roi
	roi_array[0:(n_elements(roi[*,0])-1),0:(n_elements(roi[0,*])-1),j]=roi
	xcen_array[j]=xcen & ycen_array[j]=ycen
	image=0

;	Standard grades
	if keyword_set(std) then begin
	std=where(grade eq 0)or(grade eq 2)or(grade eq 3)or(grade eq 4)or(grade eq 6)
	xs=x[std]
	ys=y[std]
	image=make_image(xs,ys)
	mkroi,image,dim+1,xcen,ycen,roi
	roi_array[0:(n_elements(roi[*,0])-1),0:(n_elements(roi[0,*])-1),j]=roi
	xcen_array[j]=xcen & ycen_array[j]=ycen
	image=0
	endif

	x=0
	y=0

endif

if(type eq 'l1')then begin
	ev=mrdfits(filename,2)
	image=fltarr(max(ev.detx)-min(ev.detx)+1,max(ev.dety)-min(ev.dety)+1)
	xarr=ev.detx-min(ev.detx) & yarr=ev.dety-min(ev.dety)
	for k=0,n_elements(yarr)-1 do begin
		image(xarr(k),yarr(k))=image(xarr(k),yarr(k)) + 1.
	endfor
	mkroi,image,dim+1,xcen,ycen,roi
	roi_array[0:(n_elements(roi[*,0])-1),0:(n_elements(roi[0,*])-1),j]=roi
	xcen_array[j]=xcen & ycen_array[j]=ycen

;	Standard grades
	if keyword_set(std) then begin
	image=fltarr(max(ev.detx)-min(ev.detx)+1,max(ev.dety)-min(ev.dety)+1)
	vec=where(ev.grade eq 0 or ev.grade eq 2 or ev.grade eq 3 $
	or ev.grade eq 4 or ev.grade eq 6)
	pis=float(ev.pi)
	pivec=pis(vec)
	for k=0,n_elements(pivec)-1 do $
	image(xarr(vec(k)),yarr(vec(k)))=image(xarr(vec(k)),yarr(vec(k))) + 1.
	mkroi,image,dim+1,xcen,ycen,roi
	roi_array[0:(n_elements(roi[*,0])-1),0:(n_elements(roi[0,*])-1),j]=roi
	xcen_array[j]=xcen & ycen_array[j]=ycen
	endif
endif

endfor
;****************************************************************
; Make postscript plots
;****************************************************************
openw,1,outfile+'_all.dat'

eer,roi_array,xcen_array,ycen_array,dim,file_array,radial_arr,ee_arr,eer_arr,wing=wing

radial,roi_array,xcen_array,ycen_array,dim,file_array,rads_arr,fav_arr,err_arr,wing=wing

pixslice,roi_array,file_array,xcen_array,ycen_array,dim,xy_arr,xz_arr,yslice_arr,zslice_arr,fwhm_arr,gauss=gauss

spec,spec_array,file_array,spec_arr,spec_range,xax_arr,gf_arr,yfit_arr,spec_out

;images2,file_array,type_array

;branch,trwid,branch_s,branch_a,branch_o

analysis_plot,file_array,type_array,roi_array,radial_arr,ee_arr,eer_arr,rads_arr,fav_arr,err_arr,xy_arr,xz_arr,yslice_arr,zslice_arr,fwhm_arr,spec_arr,xax_arr,gf_arr,yfit_arr,spec_out,spec_range



end

