pro calc_distance_counts, image, xcen, ycen, r_distance, r_counts
; Given an image array and the location of the center, this proceduce
; calculates "r_distance", an array containing the distance to each 
; pixel in the image array. It also calculates "r_counts", a second array  
; with the number of counts in each pixel. The array index in  
; r_distance and r_counts refer to the same pixel.
; VARIABLES
; image - image array
; xcen - center of image in x direction
; ycen - center of image in y  direction
; r_distance (output) - array containing the distance of each pixel in the array
;   from the center at xcen and ycen
; r_counts (output) - array containing the counts in each pixel corresponding 
;    with the r_distance array

; Author: Bob Zacher
;          Based on a routine by Brian McNamara
; Written: 1-15-98
; Property of the Smithsonian Astrophysical Observatory

array_size = size(image)
;xdim = array_size(1)
;ydim = array_size(2)

;Changed 5/12/99 to adjust max. radius
xdim=54
ydim=54

r_distance = fltarr(xdim*ydim)
r_counts = fltarr(xdim*ydim)

; For each pixel, calculate the distance from the center of the
; image

for i=0, (xdim-1) do begin

  x_distance  = float(i - xcen)
  for k=0, (ydim-1) do begin
    y_distance = float(k - ycen)
    index = k*xdim + i
    r_distance(index) = sqrt( x_distance*x_distance + y_distance*y_distance)
    ;print, r_distance(index) 
    r_counts(index) = image(i,k)

  endfor
endfor


return
end
