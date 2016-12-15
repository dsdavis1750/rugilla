pro calc_encircled_energy, image, xcen, ycen, radial_distance, cumulative_sum, encircled_energy, error,nr
; Calculates an encircled energy vector given an image and the center of
; the image 
; image - the image array
; xcen - center of image in x direction
; ycen - center of image in y  direction
; radial_distance (output) -  vector giving the positions from the center
;                             corresponding calculated encircled energy
; encircled_energy (output) - calculated  encircled energy vector
; error (output) - square root of N error bar

; Author: Bob Zacher
; Written: 1-15-98
; Property of the Smithsonian Astrophysical Observatory

calc_distance_counts, image, xcen, ycen, r_distance, r_counts

; Now calculate the cumulative sum of  the counts versus radial distance
binsize=1
cumulative_sum, r_distance, r_counts, radial_distance, cumulative_sum, unscaled_error, binsize


;--------------------------------------------------------
; changed 7/6/99 to use constant normalization radius
;  use radius of 20 pixels
;--------------------------------------------------------
;encircled_energy = cumulative_sum/ max(cumulative_sum)
; Scale the error bars properly 
;error = unscaled_error / max(cumulative_sum)

nr=20

encircled_energy = cumulative_sum / (cumulative_sum[where(radial_distance eq nr)])[0]
error = unscaled_error / (cumulative_sum[where(radial_distance eq nr)])[0]

return
end
