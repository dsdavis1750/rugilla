PRO xmm_src2reg,inst,src_file,ccf_file,plot_file,ener,frac,mllim,coor,out_file

COMMON xmmdirs,xmmdir,xmmccf

  ;open the source list file
    ;get the primary header
  phead=headfits(src_file,exten=0)
    ;get the name of the instrument
; inst=strtrim(sxpar(phead,'INSTRUME'))
    ;get the actual source list
  fxb2struct,src_file,'SRCLIST','srclist',data,ndat,head,errr
    ;do a little documentation
  dops,'p',1,10.9,10.1
  !x.range=[-2,3]
  !y.range=[0,0]
  !x.title='Log!i10!n(ML)'
  !y.title='Log!i10!n(FLux/ergs cm!e-2!n s!e-1!n)'
  plot,alog10(data.det_ml),alog10(data.flux),psym=7,symsize=0.5
  for i=0,ndat-1 do begin
    oplot,alog10(data(i).det_ml)*[1,1],alog10(data(i).flux+[-1,1]*data(i).flux_err)
  endfor
  index=findgen(51)/10.-2.0
  plc=where(data.id_band EQ 0)
  oplot,alog10(data(plc).det_ml),alog10(data(plc).flux),psym=6,color=60
; coef=poly_fit(alog10(data(plc).det_ml),alog10(data(plc).flux),1)
; oplot,index,coef(0)+coef(1)*index,color=60
  plc=where(data.id_inst EQ 0)
  oplot,alog10(data(plc).det_ml),alog10(data(plc).flux),psym=6,color=220
; coef=poly_fit(alog10(data(plc).det_ml),alog10(data(plc).flux),1)
; oplot,index,coef(0)+coef(1)*index,color=220
  oplot,alog10(mllim)*[1,1],[-20,0]
  device,/close
  spawn,'mv idl.ps '+plot_file

    ;do what selection needs to be done

  ;selecting on liklihood
  maxi=max(data.id_cluster)
  bad=intarr(maxi)
  for i=1,maxi do begin
    plc=where(data.id_cluster EQ i)
    if(plc(0) NE -1)then begin
      if(max(data(plc).det_ml) LT mllim)then bad(i-1)=1 
    endif else begin
      bad(i-1)=2
    endelse
  endfor
  bad=bad(where(bad NE 2))

  plc=where(data.id_inst EQ 0) 
  data=data(plc)
  plc=where(bad EQ 0)
  data=data(plc)
  ndat=n_elements(data)
  detx=data.detx
  dety=data.dety

  ;open the CCF file
  fxb2struct,ccf_file,'CALINDEX','calindex',indx,nndx,nhead,errr
  indx.scope=strtrim(indx.scope,2)
  indx.typeid=strtrim(indx.typeid,2)

  ;open the PSF information file
  case inst of
  "EMOS1":begin
     plc=where(indx.scope EQ 'XRT1' AND indx.typeid EQ 'XPSF')
     ccf_name=xmmccf+'XRT1_XPSF_'+$
       strtrim(string(indx(plc(0)).issue,format='(i4.4)'),2)+'.CCF'
     fxb2struct,ccf_name,'KING_PARAMS','king_params',param,npar,parhead,errr
     end
  "EMOS2":begin
     plc=where(indx.scope EQ 'XRT2' AND indx.typeid EQ 'XPSF')
     ccf_name=xmmccf+'XRT2_XPSF_'+$
       strtrim(string(indx(plc(0)).issue,format='(i4.4)'),2)+'.CCF'
     fxb2struct,ccf_name,'KING_PARAMS','king_params',param,npar,parhead,errr
     end
  "EPN":begin
     plc=where(indx.scope EQ 'XRT3' AND indx.typeid EQ 'XPSF')
     ccf_name=xmmccf+'XRT3_XPSF_'+$
       strtrim(string(indx(plc(0)).issue,format='(i4.4)'),2)+'.CCF'
     fxb2struct,ccf_name,'KING_PARAMS','king_params',param,npar,parhead,errr
     end
  endcase

    ;convert theta from rads to arcmin
  param.theta=param.theta*(180./!pi)*60.
    ;get unique ordered lists of energy and theta
  param_ener=get_unique(param.energy)
  param_thet=get_unique(param.theta)
    ;create array of rcore and alphas that correspond
  param_rcor=fltarr(n_elements(param_ener),n_elements(param_thet))
  param_alph=param_rcor
  for i=0,npar-1 do begin
    plce=where(param_ener EQ param(i).energy)
    plct=where(param_thet EQ param(i).theta)
    param_rcor(plce,plct)=param(i).params(0)
    param_alph(plce,plct)=param(i).params(1)
  endfor

  ;calculate the psfradii
  thet=xmm_off_axis_ang(detx,dety,inst)
  r_core=fltarr(ndat)
  alpha=fltarr(ndat)
  badang=intarr(ndat)
    ;find the values of r_core and alpha by interpolating from arrays
  lowe=where(param_ener LE ener)
  hige=where(param_ener GT ener)
  if(lowe(0) EQ -1 OR hige(0) EQ -1)then begin
    print,'src2reg:energy not in allowed range'
    return
  endif else begin
    lowe=max(lowe)
    hige=min(hige)
    dele=param_ener(hige)-param_ener(lowe)
    loww=(ener-param_ener(lowe))/dele
    thet_rcor=reform(param_rcor(lowe,*)*(1.-loww)+param_rcor(hige,*)*loww)
    thet_alph=reform(param_alph(lowe,*)*(1.-loww)+param_alph(hige,*)*loww)
  endelse

  for i=0,ndat-1 do begin
    lowt=where(param_thet LE thet(i))
    higt=where(param_thet GT thet(i)) 
    if(lowt(0) EQ -1 OR higt(0) EQ -1)then begin
      print,'src2reg:off-axis angle not in allowed range:',i
      badang(i)=1
    endif else begin
      lowt=max(lowt)
      higt=min(higt)
      delt=param_thet(higt)-param_thet(lowt)
      loww=(thet(i)-param_thet(lowt))/delt
      r_core(i)=thet_rcor(lowt)*(1.-loww)+thet_rcor(higt)*loww
      alpha(i) =thet_alph(lowt)*(1.-loww)+thet_alph(higt)*loww
    endelse
  endfor
  psfrad=xmm_eer(frac,r_core,alpha)
  nthe=n_elements(param_thet)
  badsub=xmm_eer(frac,thet_rcor(nthe-1),thet_alph(nthe-1))
  plc=where(badang EQ 1)
  if(plc(0) NE -1)then psfrad(plc)=badsub

  ;create the region file
  openw,olun,out_file,/get_lun
    ;print the header
  printf,olun,'# Region file format: DS9 version 3.0'
  printf,olun,'global color=green font="helvetica 10 normal" select=1 edit=1 move=1 delete=1 include=1 fixed=0 source'

  case coor of
  "sky":begin
    for i=0,ndat-1 do begin
      printf,olun,'fk5;circle('+strtrim(string(data(i).ra),2)+','+$
        strtrim(string(data(i).dec),2)+','+strtrim(string(psfrad(i)),2)+'")'
    endfor
    end
  "det":begin
    for i=0,ndat-1 do begin
      printf,olun,'physical;circle('+strtrim(string(data(i).detx),2)+','+$
        strtrim(string(data(i).dety),2)+','+strtrim(string(psfrad(i)/0.05),2)+')'
    endfor
    end
  endcase
  free_lun,olun

  openw,olun,'sas.txt',/get_lun
  case inst of
  'EMOS1':begin
    comm_str='((DETX,DETY) IN CIRCLE(100.00,-200.0,16900.))&&!('
  end
  'EMOS2':begin
    comm_str='((DETX,DETY) IN CIRCLE(0.0,-200.0,16900.))&&!('
  end
  'EPN': begin
    comm_str='!('
  end
  endcase
  for i=0,ndat-2 do begin
    comm_str=comm_str+'(DETX,DETY) IN CIRCLE('+$
      strtrim(string(data(i).detx),2)+','+$
      strtrim(string(data(i).dety),2)+','+$
      strtrim(string(psfrad(i)/0.05),2)+')||'
  endfor 
  comm_str=comm_str+'(DETX,DETY) IN CIRCLE('+$
    strtrim(string(data(i).detx),2)+','+$
    strtrim(string(data(i).dety),2)+','+$
    strtrim(string(psfrad(i)/0.05),2)+'))'
  printf,olun,comm_str
  free_lun,olun

  end

FUNCTION xmm_rad,detx,dety,inst

  case inst of
  "EMOS1":begin
    radi=sqrt((detx-100)^2+(dety+200)^2)
  end
  "EMOS2":begin
    radi=sqrt((detx-500)^2+(dety+1250)^2)
  end
  "EPN":begin
    radi=sqrt((detx-1240)^2+(dety-400)^2)
  end
  endcase

  radi=radi*0.05/60. ;answer returned in arcmin
  return,radi
  end

FUNCTION xmm_eer,frac,r_core,alpha

  rc=r_core
  al=alpha

  r=rc*sqrt(((1.-frac*(1.-((1.+((5.*60./rc)^2))^(1.-al))))^(1./(1.-al)))-1)
  return,r
  end
