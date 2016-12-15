PRO fxb2struct,file,exten,sname,data,n_entries,header,errr

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

@fxbintable

  if(n_params(0) LT 7)then begin
    print,"fxb2struct,file,exten,sname,data,n_entries,header,err"
    print,"Purpose:to read a binary FITS table into a structure"
    print,"Inputs:"
    print,"  file:the binary FITS table file"
    print,"  extension:the name or number of the extension"
    print,"  sname:the name of the structure to be created"
    print,"Outputs:"
    print,"  data:the output structure"
    print,"  n_entries:number of rows"
    print,"  header:the fits table header"
    print,"  err:error flag"
    return
  endif

  errr=0
; fxbcstruct,file,exten,sname,data,err
; if(err EQ 1)then return

  fxbopen,unit,file,exten,header
  indx=fxbfindlun(unit)
  n_entries=naxis2(indx)
  if(n_entries EQ 0)then begin
    print,'fxb2struct: no data in file'
    fxbclose,unit
    errr=1
    return
  endif
  print,'Number of table rows: ',n_entries
  n_cols=tfields(indx)

;setup for the structure

  if(n_elem(0,indx) EQ 1)then begin
    formstrg=[format(0,indx)]
  endif else begin
    formstrg=[formstrg,format(0,indx)+'('+strtrim(string(n_elem(0),indx),2)+')'] 
  endelse

  if(n_cols GT 1)then begin
    for i=2,tfields(indx) do begin
      if(n_elem(i-1,indx) EQ 1)then begin
        formstrg=[formstrg,format(i-1,indx)]
      endif else begin
        formstrg=[formstrg,format(i-1,indx)+'('+strtrim(string(n_elem(i-1,indx)),2)+')']
      endelse 
    endfor
  endif

  for i=0,tfields(indx)-1 do begin
    temp_string=ttype(i,indx)
    pos=strpos(temp_string,'-',0)
    while(pos NE -1) do begin
      temp_string=strmid(temp_string,0,pos)+'_'+strmid(temp_string,pos+1,100)
      pos=strpos(temp_string,'-',pos)
    endwhile
    ttype(i,indx)=temp_string
  endfor

  for i=0,tfields(indx)-1 do begin
    temp_string=ttype(i,indx)
    pos=strpos(temp_string,'(',0)
    while(pos NE -1) do begin
      temp_string=strmid(temp_string,0,pos)+'_'+strmid(temp_string,pos+1,100)
      pos=strpos(temp_string,'(',pos)
    endwhile
    ttype(i,indx)=temp_string
  endfor

  for i=0,tfields(indx)-1 do begin
    temp_string=ttype(i,indx)
    pos=strpos(temp_string,')',0)
    while(pos NE -1) do begin
      temp_string=strmid(temp_string,0,pos)+'_'+strmid(temp_string,pos+1,100)
      pos=strpos(temp_string,')',pos)
    endwhile
    ttype(i,indx)=temp_string
  endfor

; print,ttype

  kreate_struct,data,sname,ttype(0:n_cols-1,indx),formstrg
  data=replicate(data,n_entries)

;actually reading the structure

  if(n_entries GT 1)then begin
    for i=1,n_cols do begin
      fxbread,unit,temp,i
      data(*).(i-1)=reform(temp)
    endfor
  endif else begin
    for i=1,n_cols do begin
      fxbread,unit,temp,i
      data(*).(i-1)=temp
    endfor
  endelse

  fxbclose,unit

  end
