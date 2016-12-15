pro rad_plot,filename,energy,theta

im = readfits(filename,hdr)

s = size(im)

x = s(1)/2. 
y = s(2)/2. 

cntrd,im,x,y,xc,yc,3

print,"xc: ",xc," yc: ",yc

rads = fltarr(1,120)

radial_sub,im,x_cen,y_cen,dim,rads,fav,x,line,coeffs,slope_err, $
           minrange,maxrange,error
help,fav,rads,error
!xtitle='Radius (")'
!ytitle='Counts/pixel'
plot,rads*0.5,fav

acisi_psf_interpol,energy,theta,r,psf

s = size(fav)
eer = fltarr(s(1))
ee_error = eer

eer(0) = fav(0)
ee_error(0) = 0.

rads_nele = n_elements(rads)

for i = 1,rads_nele-1 do begin
eer(i) = eer(i-1) + fav(i)
if((eer(i-1)+fav(i)) gt 0) then begin
ee_error(i) = eer(i)*((ee_error(i-1)+error(i))/(eer(i-1)+fav(i)))
endif else begin
ee_error(i)=0.
endelse
endfor
psf = psf/max(psf)

eer = eer/eer(rads_nele-1)
;eer = eer/max(eer)

;psf = psf*(eer[1]/psf[1])

print,'max eer, max pdf',max(eer),max(psf)

ee_error=ee_error/total(eer)
;ee_error=ee_error - ee_error
!p.title=filename
!p.title=''
;plot,rads*0.5,eer,xrange=[0.1,10]
!x.range=[0.01,5]
!y.range=[0.1,1.]
ploterror,rads*0.5,eer,ee_error
!linetype=2
oplot,r,psf
!linetype=0
set_plot,'ps'
ploterror,rads*0.5,eer,ee_error
!linetype=2
oplot,r,psf
!linetype=0
psclose
!x.range=0
!y.range=0
set_plot,'X'
!p.title=' '
;stop
return
end
