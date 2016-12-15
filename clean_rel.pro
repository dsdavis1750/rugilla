PRO clean_rel,lc_file,lc_file2,extn,gt_file,smoo,plim,delt,slim,tlim,rl,time,mask,pflg

  if(n_params(0) LT 7)then begin
    print,"clean_lc,lc_file,extn,gt_file,smoo,plim,delt,slim,tlim,rl,time,mask,pflg"
    print,"Purpose:to clean a light curve"
    print,"        by creating a histogram of rate values"
    print,"        from the light curve, finding the most likely value"
    print,"        assuming that to be similar to the mean of the"
    print,"        quiescent rate, then fitting a Gaussian to a small"
    print,"        window around that value in the histogram to"
    print,"        determine the true mean and dispersion of the"
    print,"        quiescent background rate. The program then excludes"
    print,"        time intervals with rate higher than a multiple of"
    print,"        of the dispersion above the mean quiescent background"
    print,"        and excludes good regions shorter than some limit"
    print,"Inputs:"
    print,"  lc_file:name of lightcurve FITS file"
    print,"  extn:name of lightcurve extension"
    print,"  gt_file:name of output gti file"
    print,"  smoo:number of bins by which to smooth"
    print,"    should be odd"
    print,"  plim:histogram limits [min,max,step], [0.1,5.0,1.0] is good"
    print,"    plim(0)>0 to avoid problems with data drop-outs"
    print,"    plim(1)~5 is usual, program will adjust if it is too small"
    print,"    the histogram is binned by 1./smoo"
    print,"    plim(2) is a integer by which to scale histogram bin size"
    print,"  delt:+/- number of bins of window"
    print,"        around histogram peak to be fit, 10 is good for MOS"
    print,"  slim:highest sigma to be included, 2.5 is good for MOS"
    print,"  tlim:[shortest allowed good time period,"
    print,"       longest allowed bad period], [60,1] is good for MOS"
    print,"Outputs"
    print,"  rl:Gaussian fit parameters"
    print,"  time:list of times"
    print,"  mask:0 for exclude, 1 for include"
;   print,"  expr:time inclusion expression"
    print,"  pflg:1 to do plots"
    return
  endif

COMMON gauss_com,gauss_ind,gauss_dep

  plim(2)=fix(plim(2))
  fits_read,lc_file,junk,thead
  !mtitle='!17'+strtrim(string(sxpar(thead,'OBS_ID')),2)
  read_lc,lc_file,extn,time,rate
  read_lc,lc_file2,extn,tim2,rat2

  if(smoo GT 1)then begin
    smoo_rate=smooth(rate,smoo) 
    smoo_rat2=smooth(rat2,smoo)
  endif else begin
    smoo_rate=rate
    smoo_rat2=rat2
  endelse

; frac=float(n_elements(where(smoo_rate LT plim(1))))/float(n_elements(smoo_rate))
; iter=1
; while(frac LT 0.75)do begin
;   plim(1)=plim(1)*2.
;   plim(2)=plim(2)*1.
;   iter=iter+1
;   frac=n_elements(where(smoo_rate LT plim(1)))/n_elements(smoo_rate)
;   if(iter GT 4)then frac=1.0
; endwhile
; print,frac,iter,plim

  hist_cnts=histogram(smoo_rate,min=plim(0),max=plim(1),binsize=plim(2)/float(smoo))
  hist_rate=(findgen(n_elements(hist_cnts)))*plim(2)/float(smoo)+plim(0)
  maxp=where(hist_cnts EQ max(hist_cnts))
  maxp=maxp(0)
  maxr=hist_rate(maxp)

; if(smoo GT 1)then smoo_rate=smooth(rate,smoo) else smoo_rate=rate
  smoo_time=time-min(time)
  smoo_tim2=tim2-min(time)

  low=max([0,maxp-delt])
  hig=min([maxp+delt,n_elements(hist_rate)-1])

  gauss_ind=hist_rate(low:hig)
  gauss_dep=hist_cnts(low:hig)

  parms=[hist_cnts(maxp),hist_rate(maxp),0.1]
  dparm=parms/10.
  mask=[1,1,1]
  kpamoeba,'fit_gauss', parms,dparm,mask,0.01,100,rl,rh,it,cl,ch
  print,'1st ',rl,[low,hig]*plim(2)/float(smoo)
; plot,gauss_ind,gauss_dep
; stop
  ;end of the first iteration

  maxp=where(abs(hist_rate-rl(1)) EQ min(abs(hist_rate-rl(1))))
  maxp=maxp(0)
  low=max([0,maxp-delt])
  hig=min([maxp+delt,n_elements(hist_rate)-1])

  gauss_ind=hist_rate(low:hig)
  gauss_dep=hist_cnts(low:hig)
 
  parms=[hist_cnts(maxp),hist_rate(maxp),0.1]
  dparm=parms/10.
  mask=[1,1,1]
  kpamoeba,'fit_gauss', parms,dparm,mask,0.01,100,rl,rh,it,cl,ch
  print,'2nd ',rl,[low,hig]*plim(2)/float(smoo)
; plot,gauss_ind,gauss_dep
; stop
  ;end of the second iteration

  maxp=where(abs(hist_rate-rl(1)) EQ min(abs(hist_rate-rl(1))))
  maxp=maxp(0)
  low=max([0,maxp-delt])
  hig=min([maxp+delt,n_elements(hist_rate)-1])
 
  gauss_ind=hist_rate(low:hig)
  gauss_dep=hist_cnts(low:hig)
 
  parms=[hist_cnts(maxp),hist_rate(maxp),0.1]
  dparm=parms/10.
  mask=[1,1,1]
  kpamoeba,'fit_gauss', parms,dparm,mask,0.01,100,rl,rh,it,cl,ch
  print,'3rd ',rl,[low,hig]*plim(2)/float(smoo)
  ;end of the third iteration
  if(pflg)then begin
    !p.multi=[0,1,3]
    !x.omargin=[10,3]
    !y.omargin=[3,1.5]
    !x.margin=[0,0]
    !y.margin=[2,1]
    !x.range=[0,plim(1)]
    !y.range=0
    !y.title='Number'
    !x.title='Rate (cnts/s)'
    plot,hist_rate,hist_cnts,psym=10
;   oplot,hist_rate(maxp-delt)*[1,1],!y.crange
;   oplot,hist_rate(maxp+delt)*[1,1],!y.crange
    oplot,hist_rate(low)*[1,1],!y.crange
    oplot,hist_rate(hig)*[1,1],!y.crange
    oplot,hist_rate,rl(0)*exp(-(hist_rate-rl(1))^2/2./rl(2)^2),color=60

    !x.range=0
    !y.range=[0,plim(1)]
    !mtitle=''
    !y.title='Rate (cnts/s)'
    !x.title='Time -'+strtrim(string(min(time)),2)
    plot,smoo_time,smoo_rate,psym=10
    oplot,smoo_time,smoo_time*0+rl(1),color=60
    oplot,smoo_time,smoo_time*0+rl(1)+rl(2),color=220
    oplot,smoo_time,smoo_time*0+rl(1)+rl(2)*2,color=220
    oplot,smoo_time,smoo_time*0+rl(1)+rl(2)*3,color=220
  endif

;create the selection mask

  mask=smoo_time*0+1
  mask(where(smoo_rate GT rl(1)+slim*rl(2)))=0 ;simple 

;allow very short bad intervals?

  if(tlim(1) GT 0)then begin
  print,'Allowing some bad intervals'
  initpix=0
  finlpix=0
  for i=long(0),n_elements(mask)-1 do begin
    if(initpix EQ 0 and finlpix EQ 0)then begin
    ;not currently in a run
      if(mask(i) EQ 1)then begin
        finlpix=0
      endif else begin
        initpix=i
        finlpix=i
      endelse
    endif else begin
    ;currently in a run
      if(mask(i) EQ 1)then begin
        ;finished the run
        if(finlpix-initpix LT tlim(1))then begin
          mask(initpix:finlpix)=1
        endif
        initpix=0
        finlpix=0
      endif else begin
        finlpix=i
      endelse
    endelse
  endfor

  endif

;remove very short good intervals

  initpix=0
  finlpix=0
  for i=long(0),n_elements(mask)-1 do begin
    if(initpix EQ 0 and finlpix EQ 0)then begin
    ;not currently in a run
      if(mask(i) EQ 0)then begin
        finlpix=0
      endif else begin
        initpix=i
        finlpix=i
      endelse
    endif else begin
    ;currently in a run
      if(mask(i) EQ 0)then begin
        ;finished the run
        if(finlpix-initpix LT tlim(0))then begin
          mask(initpix:finlpix)=0
        endif
        initpix=0
        finlpix=0
      endif else begin
        finlpix=i
      endelse
    endelse
  endfor

  if(pflg)then begin
    oplot,smoo_time,smoo_rate*mask,psym=10,color=120
    oplot,smoo_time,smoo_time*0+rl(1),color=60
    oplot,smoo_time,smoo_time*0+rl(1)+rl(2),color=220
    oplot,smoo_time,smoo_time*0+rl(1)+rl(2)*2,color=220
    oplot,smoo_time,smoo_time*0+rl(1)+rl(2)*3,color=220

    xyouts,1000,0.8*plim(1),strtrim(string(rl(1),format='(f6.2)'),2)
    !y.title='!17Corner Rate'
    !y.range=[0,0.5]
    plot,smoo_tim2,smoo_rat2,psym=10
    oplot,smoo_tim2,smoo_rat2*mask,psym=10,color=120
    cr=total(rat2*mask)/total(mask)
    oplot,smoo_tim2,smoo_tim2*0+cr,color=220
    xyouts,1000,0.4,strtrim(string(cr,format='(f5.3)'),2)
  endif

;create the time selection list and command file

  print,'number of good intervals: ',total(mask)

  incr=5
  out_file='junk0.fits'
  cmd_file='filter_curve.csh'
  init_str="tabgtigen table="+lc_file+":RATE gtiset="+out_file
  init_str=init_str+" timecolumn=TIME expression='"

  openw,olun,cmd_file,/get_lun
 
  intr=1
  init=0
  segm=0
  rem_flg=0
  time_list='('
  out_list=''
  if(mask(0) EQ 1)then inittim=time(0) else inittim=0
  finltim=inittim
  for i=long(1),n_elements(mask)-1 do begin
    if(mask(i) EQ 1)then begin
      if(mask(i-1) EQ 1)then begin ;continue an interval
        finltim=time(i)
      endif else begin ;begin an interval
        inittim=time(i)
        finltim=time(i)
      endelse
    endif else begin
      if(mask(i-1) EQ 1)then begin ;end an interval
        if(init) then begin
          print,inittim-min(time),finltim-min(time)
          time_list=time_list+'||(TIME in ['+strtrim(string(inittim),2)+$
              ':'+strtrim(string(finltim),2)+'])'
          print,intr
          if(intr mod incr EQ 0)then begin
            init=0
            time_list=time_list+')'
            printf,olun,init_str+time_list+"'"
            out_list=out_list+' '+out_file
            time_list='('
            segm=segm+1
            out_file='junk'+strtrim(string(segm),2)+'.fits'
            init_str="tabgtigen table="+lc_file+":RATE gtiset="+out_file
            init_str=init_str+" timecolumn=TIME expression='"
          endif else begin
            init=1
          endelse
        endif else begin
          print,inittim-min(time),finltim-min(time)
          time_list=time_list+'(TIME in ['+strtrim(string(inittim),2)+$
              ':'+strtrim(string(finltim),2)+'])'
          init=1
        endelse
        if(pflg)then oplot,[inittim,finltim]-min(time),[0,0],color=20,thick=3
        intr=intr+1
      endif
    endelse
  endfor
  if(time_list NE '(')then begin
    time_list=time_list+')'
    printf,olun,init_str+time_list+"'"
    out_list=out_list+' '+out_file
  endif
  if(segm GT 0)then begin
    printf,olun,"gtimerge tables='"+out_list+"' withgtitable=yes gtitable="+gt_file+" mergemode=or" 
    printf,olun,"rm junk*.fits"
  endif else begin
    if(time_list NE '(')then begin
      printf,olun,"mv junk0.fits "+gt_file
    endif else begin
      ;serious trouble
      rem_flg=1
    endelse
  endelse
  free_lun,olun
  if(rem_flg)then spawn,'rm filter_curve.csh'

  print,'number of good intervals: ',total(mask)
  print,'of: ',n_elements(mask)
  end
