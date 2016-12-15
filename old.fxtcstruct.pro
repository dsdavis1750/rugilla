PRO fxtcstruct,file,exten,sname,struct,err

;This routine reads a FITS table and creates an IDL structure.
;The column names in the table are the structure tags
;The structure will include values from the first row of the table

;inputs:
;  file:name of the FITS file
;  exten:number of the extension containing the binary table
;  sname:name of the structure to create
;outputs:
;  struct:the structure
;  err:an error flag

  err=0
  openw,olun,'temp_fun.pro',/get_lun
  printf,olun,'FUNCTION temp_fun,struct'
  printf,olun,' '
  stst="  struct=kreate_struct(name='"+sname+"',"
  nmst="["
  vlst=" "
  fits_read,file,tab,htab,exten=exten
  ftinfo,htab,hstruct

  nfields=n_elements(hstruct.ttype)
  for i=1,nfields do begin
    nmst=nmst+"'"+strtrim(hstruct.ttype(i-1),2)+"',"
    data=ftget(hstruct,tab,i)
    data=max(data)
    tform=strmid(hstruct.tform(i-1),0,1)
    if(tform EQ 'A')then vlst=vlst+"'"
    vlst=vlst+strtrim(string(data),2)
    if(tform EQ 'A')then vlst=vlst+"'"
    vlst=vlst+','
  endfor
  len=strlen(vlst)
  vlst=strmid(vlst,0,len-1)
  len=strlen(nmst)
  nmst=strmid(nmst,0,len-1)

  stst=stst+nmst+"],"+vlst+")"

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
