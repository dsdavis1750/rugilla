PRO read_pha,spec_file,chan,cnts,stat,qual,grup,expo,back,area,err,err_msg

;stat added 24 June 2002

  if(n_params(0) LT 9)then begin
    print,"read_pha,spec_file,chan,cnts,stat,qual,grup,expo,back,area,err"
    print,"Purpose:reads a OGIP style PHA/PI file"
    print,"Inputs:"
    print,"  spec_file:name of the spectrum file"
    print,"Outputs:"
    print,"  chan:vector of channel numbers"
    print,"  cnts:vector of counts"
    print,"  stat:vector of statistical errors"
    print,"  qual:vector of quality flags"
    print,"  grup:vector of group flags"
    print,"  expo:exposure time (scalor)"
    print,"  back:backscal parameter (scalor)"
    print,"  area:areascal parameter (scalor)"
    print,"  err:error flag"
    return
  endif

  err=0
  err_msg=' '

  needed=['CHANNEL','COUNTS','QUALITY','GROUPING']
  n_needed=n_elements(needed)

  fxbopen,unit,spec_file,'spectrum',h ;open the corrext binary table'
  fxbfind,unit,'ttype',cols,vals,nfnd ;get the names of the columns

  vals=strupcase(strtrim(vals,2))
 
  plc=where(vals EQ 'CHANNEL')
  if(plc(0) EQ -1)then begin
    err=1
    err_msg='No channel column'
    return
  endif else begin
    fxbread,unit,chan,plc(0)+1
  endelse

  rflg=0
  plc=where(vals EQ 'COUNTS')
  if(plc(0) EQ -1)then begin
    plc=where(vals EQ 'RATE')
    if(plc(0) EQ -1)then begin
      err=1
      err_msg='No counts column'
      return
    endif else begin
      rflg=1
      fxbread,unit,cnts,plc(0)+1
    endelse
  endif else begin
    fxbread,unit,cnts,plc(0)+1
  endelse

  plc=where(vals EQ 'STAT_ERR')
  if(plc(0) EQ -1)then begin
    if(sxpar(h,'POISSERR') EQ 1)then begin
      stat=sqrt(abs(cnts))
    endif
  endif else begin
    fxbread,unit,stat,plc(0)+1
  endelse

  plc=where(vals EQ 'QUALITY')
  if(plc(0) EQ -1)then begin
    qual=chan*0
  endif else begin
    fxbread,unit,qual,plc(0)+1
  endelse

  plc=where(vals EQ 'GROUPING')
  if(plc(0) EQ -1)then begin
    grup=chan*0+1
  endif else begin
    fxbread,unit,grup,plc(0)+1
  endelse

  fxbclose,unit

  expo=sxpar(h,'EXPOSURE')
  if(!err EQ -1)then begin
    err=1
    err_msg='No exposure keyword found'
    return
  endif 
  area=sxpar(h,'AREASCAL')
  if(!err EQ -1)then begin
    err=1
    err_msg='No areascal keyword found'
    return
  endif 
  back=sxpar(h,'BACKSCAL')
  if(!err EQ -1)then begin
    err=1
    err_msg='No backscal keyword found'
    return
  endif 

  back=float(back)
  area=float(area)
  end
