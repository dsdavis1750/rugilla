pro cumulative_sum, x_val_array, counts_array, bin_array, cumulative_sum, error, binsize
; Produces a cumulative sum given an array containing a list of distances 
; and a corresponding array containing counts
; x_val_array - array with list of distances
; counts_array - contains counts corresponding to distances in x_val_array
; bin_array - (output) the position array for the output cumulative sum
; cumulative_sum - (output) the cumulative sum corresponding to the positions
;                  in bin_array
; error - (output) the square root of N error on cumulative_sum
; binsize - the desired binning size for the output arrays
; Author: Bob Zacher
; 1-15-98
; Property of the Smithsonian Astrophysical Observatory


; find the maximum  x value
x_max = max(x_val_array)
out_array_size = fix(x_max+1)/ binsize
bin_array = findgen( out_array_size ) ; distance from center
binned_counts = fltarr( out_array_size) ; total counts in an annulus of 
                                        ; thickness binsize
cumulative_sum = fltarr( out_array_size-1)     ; cumulative sum up to a given radius
error = fltarr( out_array_size-1 )
num_events = n_elements( x_val_array )

for i=0, (out_array_size-2) do begin

  indices = where( (x_val_array ge bin_array(i))  $
       and (x_val_array lt bin_array(i+1)))

  if ( indices(0) ne -1) then begin  
    result = counts_array(indices)
    binned_counts(i) = total( result )

  endif else begin
    result = 0
    binned_counts(i) = 0
  endelse

  ; Calculate the cumulative sum

  if ( i eq 0 ) then begin
    cumulative_sum(i)  = 0
  endif else begin
    cumulative_sum(i)  = cumulative_sum(i-1) + binned_counts(i-1)
    ; Now calculate the error. Take the square root of the cumulative counts 
    error(i) = sqrt( cumulative_sum(i) ) 
  endelse


endfor


return
end
