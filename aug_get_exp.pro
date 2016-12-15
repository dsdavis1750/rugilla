PRO aug_get_exp,db_file,c_val,c_dst,c_tim,p_flg,aug,aug_exp

COMMON spec_data,emin,emax,temp_spec,temp_uncr,temp_resu

  if(n_params(0) LT 7)then begin
    print,"get_aug,db_file,r_lim,h_lim,d_lim,c_val,c_dst,c_tim,p_flg,aug,aug_exp"
    print,"Purpose:to create a mean corner spectrum"
    print,"  Given a rate, hardness, and date,"
    print,"  finds the closest c_tim seconds of exposure"
    print,"  in that 3-D space"
    print,"Inputs:"
    print,"  db_file:database file of corner spectra"
    print,"  c_val:central [rate,hardness,date]"
    print,"  c_dst:distance scaling in those directions"
    print,"    as factor of mean uncertainty in those directions"
;   print,"    suggest:[0.00323,0.358,1.0]"
    print,"  c_tim:how many seconds of exposure desired for spectrum"
    print,"  p_flg:1 if informative plot to be made"
    print,"Outputs:"
    print,"  aug:the augmented spectrum"
    print,"  aug_exp:the amount of exposure time in the augmented spectrum"
    return
  endif

  fxbopen,unit,db_file,'QUIES',qhead
  fxbread,unit,obid,'OBSID'
  fxbread,unit,name,'NAME'
  fxbread,unit,objt,'OBJECT'
  fxbread,unit,date,'DATE'
  fxbread,unit,revl,'REVOLUT'
  fxbread,unit,filt,'FILTER'
  fxbread,unit,time,'ONTIME'
  fxbread,unit,expo,'EXPO'
  fxbread,unit,rate,'RATE'
  fxbread,unit,urate,'RATE_ERR'
  fxbread,unit,hard,'HARD'
  fxbread,unit,uhard,'HARD_ERR'
  fxbread,unit,spec,'SPECTRUM'
  fxbclose,unit

;do some cleaning

  plc=where(finite(rate) and finite(hard) and expo GT 1000.)
  obid=obid(plc)
  name=name(plc)
  objt=objt(plc)
  date=date(plc)
  revl=revl(plc)
  filt=filt(plc)
  time=time(plc)
  expo=expo(plc)
  rate=rate(plc)
  urate=urate(plc)
  hard=hard(plc)
  uhard=uhard(plc)
  spec=spec(*,plc)

  nobs=n_elements(plc)
  energ=findgen(800)*15./1000.

;get some temporary values

  murate=total(urate)/n_elements(urate)
  muhard=total(uhard)/n_elements(uhard)
  print,'Rate: ',c_val(0),' Rate error: ',murate
  print,'Hard: ',c_val(1),' Hard error: ',muhard

;find the closest c_dst(0) secs

  if(c_tim EQ 0)then begin
    print,'Error: no exposure length provided'
    return
  endif else begin
    if(c_dst(0) EQ 0 and c_dst(1) EQ 0 AND c_dst(2) EQ 0)then begin
      print,'Error: no distance scalings'
      return
    endif else begin
      good=intarr(nobs)
      if(c_dst(0) NE 0)then begin
        rdist=(rate-c_val(0))/murate/c_dst(0) 
      endif else begin
        rdist=rate*0
      endelse

      if(c_dst(1) NE 0)then begin
        hdist=(hard-c_val(1))/muhard/c_dst(1) 
      endif else begin
        hdist=hard*0
      endelse

      if(c_dst(2) NE 0)then begin
        ddist=(revl-c_val(2))/c_dst(2) 
      endif else begin
        ddist=revl*0
      endelse

      dist=sqrt(rdist^2+hdist^2+ddist^2)
      ord=sort(dist)
      cum_exp=0
      i=0
      while(cum_exp LT c_tim and i LT nobs)do begin
        cum_exp=cum_exp+expo(ord(i))
        good(ord(i))=1
;       print,i,ord(i),expo(ord(i)),cum_exp
        i=i+1
      endwhile
    endelse
  endelse

  good=where(good NE 0)

  if(p_flg)then aug_plot,rate,urate,hard,uhard,revl,good,c_val

;actually put together the spectrum

  aug=fltarr(800)
  aug_exp=0.
  for i=0,n_elements(good)-1 do begin
    aug=aug+spec(*,good(i)) 
    aug_exp=aug_exp+expo(good(i))
  endfor

  end  
