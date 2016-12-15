PRO read_lc,file,exten,time,rate

  if(n_params(0) LT 3)then begin
    print,"read_lc,file,exten,time,rate"
    print,"Purpose:read a lightcurve file"
    print,"Inputs:"
    print,"  file:name of lightcurve file"
    print,"  exten:name of the lightcurve extension"
    print,"Outputs:"
    print,"  time:vector of times"
    print,"  rate:vector of rates"
    return
  endif

  fxbopen,unit,file,exten,head
  fxbfind,unit,'TTYPE',cols,vals,n_found
  vals=strupcase(vals)
  time_col=cols(where(strtrim(vals,2) EQ 'TIME'))
  time_col=time_col(0)
  rate_col=cols(where(strtrim(vals,2) EQ 'RATE'))
  rate_col=rate_col(0)
  fxbread,unit,time,time_col
  fxbread,unit,rate,rate_col
  fxbclose,unit

  end

