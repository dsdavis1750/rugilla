pro format,dummy

ftime = 7.428703703703703D-4 

form = '(a9,1x,e21.15)'
form1 = '(e21.15)'

;print,format=form,'MJDREFF = ',ftime

txt = 'MJDREF keyword comments'

; write a fits file

image = intarr(10,10)

MKHDR, header, image

SXADDPAR, header, 'MJDREFF', ftime, 'MJDREF keyword comments', FORMAT=form1

print,header

return
end


