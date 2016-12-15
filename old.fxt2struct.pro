PRO fxt2struct,file,exten,sname,data,n_entries,header,err

;This routine reads a FITS binary table into and IDL structure.
;The column names in the table are the structure tags

;inputs:
;  file:name of the FITS binary table
;  exten:the extension number containing the binary table
;  sname:the name of the structure
;outputs:
;  data:the structure containing the data
;  n_entries:the number of entries (rows) in the table
;  err:an error flag

  fxtcstruct,file,exten,sname,data,err
  if(err EQ 1)then return

  fits_read,file,tab,htab,exten=exten
  ftinfo,htab,hstruct

  n_entries=sxpar(htab,'NAXIS2')
  print,'Number of table rows: ',n_entries
  n_cols=n_elements(hstruct.tform)

  if(n_entries GT 0)then begin
    tems=replicate(data,n_entries)
    for i=1,n_cols do begin
      tems(*).(i-1)=ftget(hstruct,tab,i)
    endfor
    data=tems
  endif

  end
