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

; fxtcstruct,file,exten,sname,data,err
; if(err EQ 1)then return

  htab=headfits(file,exten=exten)
  ftinfo,htab,hstruct
  n_entries=sxpar(htab,'naxis2')
  if(n_entries EQ 0)then begin
    errr=1
    return
  endif
  print,'Number of table rows: ',n_entries
  n_cols=n_elements(hstruct.tbcol)

;setup for the structure

  formstrg=strarr(n_cols)
  for i=0,n_cols-1 do formstrg(i)=strmid(hstruct.tform(i),0,1)
  for i=0,n_cols-1 do begin
    plc=where(hstruct.ttype EQ hstruct.ttype(i))
    if(n_elements(plc) GT 1)then begin
      hstruct.ttype(i)=hstruct.ttype(i)+'0'
    endif
  endfor

  for i=0,n_cols-1 do begin
    temp_string=hstruct.ttype(i)
    pos=strpos(temp_string,'-',0)
    while(pos NE -1) do begin
      temp_string=strmid(temp_string,0,pos)+'_'+strmid(temp_string,pos+1,100)
      pos=strpos(temp_string,'-',pos)
    endwhile
    hstruct.ttype(i)=temp_string
  endfor

  for i=0,n_cols-1 do begin
    temp_string=hstruct.ttype(i)
    pos=strpos(temp_string,'(',0)
    while(pos NE -1) do begin
      temp_string=strmid(temp_string,0,pos)+'_'+strmid(temp_string,pos+1,100)
      pos=strpos(temp_string,'(',pos)
    endwhile
    hstruct.ttype(i)=temp_string
  endfor

  for i=0,n_cols-1 do begin
    temp_string=hstruct.ttype(i)
    pos=strpos(temp_string,')',0)
    while(pos NE -1) do begin
      temp_string=strmid(temp_string,0,pos)+'_'+strmid(temp_string,pos+1,100)
      pos=strpos(temp_string,')',pos)
    endwhile
    hstruct.ttype(i)=temp_string
  endfor
  
  kreate_struct,data,sname,hstruct.ttype,formstrg

;then read the data

  fits_read,file,tab,htab,exten=exten

  if(n_entries GT 0)then begin
    tems=replicate(data,n_entries)
    for i=1,n_cols do begin
      tems(*).(i-1)=ftget(hstruct,tab,i)
    endfor
    data=tems
  endif

  end
