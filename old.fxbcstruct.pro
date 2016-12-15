PRO fxbcstruct,file,exten,sname,struct,err

;This routine reads a binary FITS table and creates an IDL structure.
;The column names in the table are the structure tags
;The structure will include values from the first row of the table

;inputs:
;  file:name of the FITS file
;  exten:number of the extension containing the binary table
;  sname:name of the structure to create
;outputs:
;  struct:the structure
;  err:an error flag

@fxbintable

  err=0
  openw,olun,'temp_fun.pro',/get_lun
  printf,olun,'FUNCTION temp_fun,struct'
  printf,olun,' '
  stst="  struct=kreate_struct(name='"+sname+"',"
  nmst="["
  vlst=" "
  fxbopen,unit,file,exten,h
  indx=fxbfindlun(unit)
  if(naxis2(indx) EQ 0)then begin
    fxbclose,unit
    free_lun,olun
    err=1
    return
  endif

  nfields=tfields(indx)
  for i=1,nfields do begin
    nmst=nmst+"'"+strtrim(ttype(i-1,indx),2)+"',"
    fxbread,unit,data,i,1
    if(n_elem(i-1,indx) NE 1)then begin
      vlst=vlst+'['
      for j=0,n_elem(i-1,indx)-2 do begin
        if(format(i-1,indx) EQ 'A' or format(i-1,indx) EQ 'L')then vlst=vlst+"'"
        if(format(i-1,indx) EQ 'E' AND NOT finite(data(j)))then data(j)=0.0 
        vlst=vlst+strtrim(string(data(j)),2)
        if(format(i-1,indx) EQ 'A' or format(i-1,indx) EQ 'L')then vlst=vlst+"'"
        vlst=vlst+','
      endfor
      if(format(i-1,indx) EQ 'A' or format(i-1,indx) EQ 'L')then vlst=vlst+"'"
      if(format(i-1,indx) EQ 'E' AND NOT finite(data(j)))then data(j)=0.0 
      vlst=vlst+strtrim(string(data(j)),2)
      if(format(i-1,indx) EQ 'A' or format(i-1,indx) EQ 'L')then vlst=vlst+"'"
      vlst=vlst+'],'
    endif else begin
      if(format(i-1,indx) EQ 'A' or format(i-1,indx) EQ 'L')then vlst=vlst+"'"
      if(format(i-1,indx) EQ 'E')then begin
;       print,format(i-1,indx),i,indx,data
        if( NOT finite(data))then begin
          data=0.0 
        endif
      endif
      vlst=vlst+strtrim(string(data),2)
      if(format(i-1,indx) EQ 'A' or format(i-1,indx) EQ 'L')then vlst=vlst+"'"
      vlst=vlst+','
    endelse
  endfor
  len=strlen(vlst)
  vlst=strmid(vlst,0,len-1)
  len=strlen(nmst)
  nmst=strmid(nmst,0,len-1)

  stst=stst+nmst+"],"+vlst+")"

  fxbclose,unit

  result=execute(stst)

  tlen=strlen(stst)
  while(tlen GT 75)do begin
    temp=strmid(stst,0,75)
    j=74
    while(strmid(stst,j,1)NE',')do j=j-1
    printf,olun,'  '+strmid(stst,0,j+1)+'$'
    stst=strmid(stst,j+1,1000)
    tlen=strlen(stst)
  endwhile
  printf,olun,'  '+stst
  printf,olun,' '
  printf,olun,'  return,a'
  printf,olun,'  end'

  free_lun,olun

;  struct=call_function('temp_fun')
  end
