PRO fxb2estruct,file,exten,data,n_entries,header,err

;This routine reads a FITS binary table into and IDL structure.
;The column names in the table are the structure tags

;inputs:
;  file:name of the FITS binary table
;  exten:the extension number containing the binary table
;  data:the structure to which the data is to be transferred
;outputs:
;  n_entries:the number of entries (rows) in the table
;  err:an error flag

@fxbintable

  fxbopen,unit,file,exten,header
  indx=fxbfindlun(unit)
  n_entries=naxis2(indx)
  n_cols=tfields(indx)

  if(n_entries GT 0)then begin
    tems=replicate(data,n_entries)
    for i=1,n_cols do begin
      fxbread,unit,temp,i
      tems(*).(i-1)=temp
    endfor
    data=tems
  endif

  fxbclose,unit

  end
