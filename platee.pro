PRO plate3,mosnum,obj_file,obj_rmf,bck_file,lim_area,pflag,errr

;this program is somewhat intelligent, returns an error if a file is missing

COMMON xmmdirs,xmmdir

  if(n_params(0) LT 3)then begin
    print,"platee,mosnum,obj_file,obj_rmfi,bck_file,errr"
    print,"Purpose:to construct the MOS particle background"
    print,"Inputs:"
    print,"  mosnum:either 1 or 2"
    print,"  obj_file:name of object spectrum file"
    print,"  obj_rmf:RMF file for the object region"
    print,"  bkg_file:name of output object background spectrum file"
    print,"  lim_area:value of backscale below which chips are eliminated"
    print,"    300000. is a good number to use"
    print,"  pflag:1 for plotting"
    return
  endif

  print,'WARNING: if chip 3, 4 or 6 are missing spurious features can arise'
  print,'         in 2.04-2.33, 9.38-9.98, 11.03-11.82, the gold regions'

;modified on 18 Jan 04 to fix minor problems
;modified on 20 Oct 04 to allow for missing files
;modified on 20 Oct 04 to properly mask out the gold features
;modified on 19 Oct 04 to allow arbitrary spectrum length
;modified on  4 Aug 04 to allow one to discard chips with too little coverage

;the background for the region of interest =
;sum over all of the chips of (observation corner)*(FWC in region)
;                                           (FWC corner)
;the complications:
;  augmenting the corner data based on rate and hardness
;  missing chips 
;  chip 1, which has no corner data
;    solve by averaging "similar chips"
;  features in chip 2,7 corner data not present in FOV

  sls_col,0

  errr=0
  aug_time=1e6
  smo_fact=9
  aug_scal=[1,1,0]

;some setup

  case mosnum of
  1:begin
    sim_chips=[2,3,6,7]
    n_sims=4
    mosnam='1'
    end
  2:begin
    sim_chips=[2,4,6,7]
    n_sims=4
    mosnam='2'
    end
  endcase

  if(!d.name EQ 'PS')then begin
    dops,'p',1,10.9,24.1
  endif

;get information about the object spectrum

  read_pha,obj_file,chan,obj_cnts,stat,qual,grup,obj_expo,obj_back,area,err
  objc=float(obj_cnts)/obj_back/obj_expo
  n_chan=n_elements(chan)

;get the channel to energy conversion
  
  read_rmf,obj_rmf,elow,ehig,detm,chan,emin,emax
  ener=(emin+emax)/2.
  plc_al=where(ener GE 0.3 AND ener LE 10.0) ;region over which background cald
  plc_lo=where(ener GE 0.4 AND ener LE 0.8) ;a soft band
  plc_hi=where(ener GE 2.5 AND ener LE 5.0) ;a hard band
  ;places to remove the gold lines
  plc_au=where((ener GE 2.04 and ener LE 2.33) $
              OR (ener GE 9.38 and ener LE 9.98) $
              OR (ener GE 11.03 and ener LE 11.82))
  ;the place to remove the aluminum and silicon lines
  plc_alsi=where(ener GE 1.2 and ener LE 1.90)

;read the fwc spectra from the object region
;  figure out which chips are involved
;  use reg_area as a vector describing which chips are in use

  reg_spec=fltarr(n_chan,7)
  reg_uncr=fltarr(n_chan,7)
  reg_area=fltarr(7)

  reg_rati=fltarr(n_chan,7)
  reg_mask=intarr(n_chan,7)

  for i=1,7 do begin
    spec_name='mos'+mosnam+'-'+string(i,format='(i1)')+'ff.pi'
    test_file=findfile(spec_name)
    if(test_file EQ spec_name)then begin ;the file exists
      fits_read,spec_name,junk,thead
      fxbopen,unit,spec_name,'SPECTRUM',head
      fxbread,unit,spec,2
      fxbclose,unit

      spec=float(spec)
      expo=sxpar(head,'EXPOSURE')
      back=sxpar(thead,'BACKSCAL')
      if(back NE 0 and expo NE 0)then begin ;there is active area on this chip
        ;normalize
        uncr=sqrt(spec)
        reg_spec(*,i-1)=spec/back/expo
        reg_uncr(*,i-1)=uncr/back/expo
        reg_area(i-1)=back
      endif
    endif else begin
      reg_spec(*,i-1)=0.0
      reg_uncr(*,i-1)=0.0
      reg_area(i-1)=0.0
    endelse
  endfor

;remove any chip that has less area than lim_area
  for i=0,6 do begin
    print,'Chip ',strtrim(string(i),2),' backscale=',$
      strtrim(string(reg_area(i)),2)
  endfor
  plc=where(reg_area LT lim_area)
  if(plc(0) NE -1)then begin
    for i=0,n_elements(plc)-1 do begin
      reg_area(plc(i))=0
      reg_spec(*,plc(i))=0
      reg_uncr(*,plc(i))=0
    endfor
  endif

;read corner spectra from the observation

  oc_spec=fltarr(n_chan,7)
  oc_uncr=fltarr(n_chan,7)
  oc_mask=fltarr(n_chan,7)+1
  oc_mask(plc_au,1)=0
  oc_mask(plc_au,6)=0
  oc_area=fltarr(7)

  if(pflag)then begin
    !p.multi=[0,1,3]
    !x.range=0
    !y.range=0
  endif

  ;even if one of the outer chips is not being used
  ;we still need the corner data if chip 1 is on
  for i=2,7 do begin
    spec_name='mos'+mosnam+'-'+string(i,format='(i1)')+'oc.pi'
    test_file=findfile(spec_name)
    if(test_file EQ spec_name)then begin ;the file exists
      fits_read,spec_name,junk,thead
      fxbopen,unit,spec_name,'SPECTRUM',head
      fxbread,unit,spec,2
      fxbclose,unit
      spec=float(spec)
      expo=sxpar(head,'EXPOSURE')
      revl=sxpar(thead,'REVOLUT')
      back=sxpar(thead,'BACKSCAL')

      rate=total(spec(plc_al))/expo
      hard=total(spec(plc_hi))/total(spec(plc_lo))

    ;augment
      db_file=xmmdir+'qpb/mos'+mosnam+'-'+string(i,format='(i1)')+'.fits'
      aug_get_exp,db_file,[rate,hard,revl],aug_scal,aug_time,$
        pflag,aug_spec,aug_expo

      spec=spec+aug_spec
      expo=expo+aug_expo
      aug_rate=total(spec(plc_al))/expo
      spec=spec*rate/aug_rate
   
    ;normalize
      uncr=sqrt(spec)
      oc_spec(*,i-1)=spec/back/expo
      oc_uncr(*,i-1)=uncr/back/expo
      oc_area(i-1)=back 
    endif else begin
      oc_spec(*,i-1)=0.0
      oc_uncr(*,i-1)=0.0
      oc_area(i-1)=0.0
      if(reg_area(i-1) NE 0)then begin
        errr=1
        reg_spec(*,i-1)=0.0
        reg_uncr(*,i-1)=0.0
        reg_area(i-1)=0
        print,'missing file: ',spec_name,' in platee'
        ;but we won't return, we'll try to recover
      endif
    endelse
  endfor

  if(pflag and !d.name EQ 'PS')then begin
    device,/close
    spawn,'mv idl.ps augment.ps'
  endif

;for the center chip
  if(reg_area(0) NE 0)then begin ;the center chip is active
    ;create a mask for normalizing the center chips spectrum
    mk1=intarr(n_chan)

    ;determine how many similar chips have data
    n_real=n_sims
    for i=0,n_sims-1 do begin
      if(oc_area(sim_chips(i)-1) EQ 0)then begin
        n_real=n_real-1
        sim_chips(i)=0
      endif
    endfor    
    if(n_real EQ 0)then begin
      ;this is really an error from which we can't recover
      errr=2
      print,'No data for chip 1 corners'
      return
    endif else begin
      n_sims=n_real
      plc=where(sim_chips NE 0)
      sim_chips=sim_chips(plc)
    endelse
    print,'N_SIMS: ',n_sims
    print,sim_chips

    ;sum over the similar chips
    oc1=fltarr(n_chan)
    uc1=fltarr(n_chan)
    for i=0,n_sims-1 do begin
      oc1=oc1+oc_spec(*,sim_chips(i)-1)*oc_mask(*,sim_chips(i)-1)
      uc1=uc1+(oc_spec(*,sim_chips(i)-1)^2)*oc_mask(*,sim_chips(i)-1)
      mk1=mk1+oc_mask(*,sim_chips(i)-1)
    endfor
    plc=where(mk1 NE 0)
    if(plc(0) EQ -1)then begin
      errr=2
      print,'another problem'
      return
    endif
    oc1(plc)=oc1(plc)/mk1(plc)
    uc1(plc)=sqrt(uc1(plc))/mk1(plc)
    plc=where(mk1 EQ 0)
    if(plc(0) NE -1)then begin
      oc1(plc)=0.0
      uc1(plc)=0.0
    endif
    oc_spec(*,0)=oc1
    oc_uncr(*,0)=uc1
  endif

;read FWC corner spectra

  if(pflag)then !p.multi=[0,1,7]

  fc_spec=fltarr(n_chan,7)
  fc_uncr=fltarr(n_chan,7)
  fc_area=fltarr(7)

  for i=2,7 do begin
    spec_name='mos'+mosnam+'-'+string(i,format='(i1)')+'fc.pi'
    test_file=findfile(spec_name)
    if(test_file EQ spec_name)then begin ;the file exists
      fits_read,spec_name,junk,thead
      fxbopen,unit,spec_name,'SPECTRUM',head
      fxbread,unit,spec,2
      fxbclose,unit
       spec=float(spec)
      expo=sxpar(head,'EXPOSURE')
      back=sxpar(thead,'BACKSCAL')

    ;normalize
      uncr=sqrt(spec)
      fc_spec(*,i-1)=spec/back/expo
      fc_uncr(*,i-1)=sqrt(spec)/back/expo
      fc_area(i-1)=back
    endif else begin
      if(reg_area(i-1) NE 0 OR oc_area(i-1) NE 0)then begin
        errr=1
        print,'missing file: ',spec_name,' in platee'
        return
      endif
    endelse
  endfor

  ;for the center chip
  if(reg_area(0) NE 0)then begin ;newest change 15/11/04
    fc1=fltarr(n_chan)
    for i=0,n_sims-1 do begin
      fc1=fc1+fc_spec(*,sim_chips(i)-1)*oc_mask(*,sim_chips(i)-1)
    endfor
    plc=where(mk1 NE 0)
    if(plc(0) EQ -1)then begin
      errr=2
      print,'yet another problem'
      return
    endif
    fc1(plc)=fc1(plc)/mk1(plc)
    plc=where(mk1 EQ 0)
    if(plc(0) NE -1)then fc1(plc)=0.0
    fc_spec(*,0)=fc1
  endif ;newest change
  
;create the ratio of region/corner for FWC data

  for i=2,7 do begin
    if(reg_area(i-1) NE 0)then begin
      temp_top=smooth(reg_spec(*,i-1),smo_fact)
      temp_bot=smooth(fc_spec(*,i-1),smo_fact)
      plc=where(temp_bot NE 0)
      reg_mask(plc,i-1)=1.
      reg_rati(plc,i-1)=temp_top(plc)/temp_bot(plc) 
    endif
  endfor

  ;single out chip 1 for special treatment again
  if(reg_area(0) NE 0)then begin
    ;does the spectrum have gold gaps?
    if(min(oc_spec(plc_al,0)) NE 0)then begin
      temp_top=smooth(reg_spec(*,0),smo_fact)
      temp_bot=smooth(fc_spec(*,0),smo_fact)
      plc=where(temp_bot NE 0)
      reg_mask(plc,0)=1.
      reg_rati(plc,0)=temp_top(plc)/temp_bot(plc) 
    endif else begin
      errr=3
      print,'You are scrod, your chip 1 background has goldgaps'
      return
    endelse
  endif

;some documentation type stuff

  if(pflag)then begin
    if(!d.name EQ 'X')then window,0,xsize=512,ysize=1024
    !p.multi=[0,1,7]
    !x.range=[0.3,12]
    !y.range=[1e-13,2e-11]
    for i=1,7 do begin
      !mtitle='!17MOS '+mosnam+' Chip '+strtrim(string(i),2)
      fact=total(fc_spec(plc_hi,i-1))/total(oc_spec(plc_hi,i-1))
      ;first, plot the fwc spectrum for object region (black)
      plot_oo,ener(plc_al),smooth(reg_spec(plc_al,i-1),9),psym=10
      ;then plot the fwc corner spectrum (dark blue)
      oplot,ener(plc_al),smooth(fc_spec(plc_al,i-1),9),psym=10,color=60
      ;now, to give some idea of what the noise is in this procedure,
      ;fwc corner*ratio = (ideally) fwc fov (red)
      oplot,ener(plc_al),fc_spec(plc_al,i-1)*reg_rati(plc_al,i-1),$
        psym=10,color=220
      ;one hopes that the following will be as smooth as possible (green)
      oplot,ener(plc_al),ener(plc_al)*0+2e-13
      oplot,ener(plc_al),reg_rati(plc_al,i-1)*2e-12,psym=10,color=120
      ;the following will not be equal to the red line because the corner spectra
      ;are derived from different time periods, and thus have different shapes.
      ;this is the background spectrum to be used for this chip (orange)
      oplot,ener(plc_al),smooth(oc_spec(plc_al,i-1),9)*reg_rati(plc_al,i-1)*$
        fact,psym=10,color=180
      oplot,ener(plc_al),oc_spec(plc_al,i-1)*reg_rati(plc_al,i-1)*fact,$
        psym=10,color=180
    endfor
    if(!d.name EQ 'PS')then begin
      device,/close
      spawn,'mv idl.ps chips.ps'
    endif
  endif

  print,fc_area
  print,reg_area
  print,total(fc_area),total(reg_area),$
    total(fc_area)/(total(fc_area)+total(reg_area))

;create the masks to remove strong instrumental lines

  ;remove the Al-Si line from all spectra
  if(plc_alsi(0) NE -1)then oc_mask(plc_alsi,*)=0
  ;remove the Au lines from chips 2&7
    ;already done above

;sum the background over all active chips

  spec=fltarr(n_chan)
  uncr=fltarr(n_chan)
  mask=fltarr(n_chan)
  for i=1,7 do begin
    if(reg_area(i-1) NE 0)then begin
      plc=where(reg_mask(*,i-1) NE 0 and oc_mask(*,i-1) NE 0)
      spec(plc)=spec(plc)+oc_spec(plc,i-1)*reg_rati(plc,i-1)*reg_area(i-1)
      uncr(plc)=uncr(plc)+(oc_uncr(plc,i-1)*reg_rati(plc,i-1)*reg_area(i-1))^2
      mask(plc)=mask(plc)+reg_area(i-1)
    endif
  endfor  

  plc=where(mask NE 0)
  if(plc(0) NE -1)then begin
    spec(plc)=spec(plc)/mask(plc)
    uncr(plc)=sqrt(uncr(plc))/mask(plc)
  endif

;bridge the Al-Si gap
  ;set a number of line-free windows
  ;extract the data from those windows
  ;fit a simple polynomial
  ;use to interpolate over the Al-Si gap

  ;the following two lines need to be converted to energies
  mine=[40,127,256,370,405,440,475,505,545,585,665,787]
  maxe=[79,255,350,385,415,460,490,525,560,625,735,799]
  nels=n_elements(maxe)
  tempe=[ener(mine(0):maxe(0))]
  temps=[spec(mine(0):maxe(0))]
  tempu=[uncr(mine(0):maxe(0))]
  for i=1,nels-1 do begin
    tempe=[tempe,ener(mine(i):maxe(i))]
    temps=[temps,spec(mine(i):maxe(i))]
    tempu=[tempu,uncr(mine(i):maxe(i))]
  endfor
  plc=where(tempe GT 0 and temps GT 0)
  tempe=alog(tempe(plc))
  temps=alog(temps(plc))
  nord=5
  coef=poly_fit(tempe,temps,nord)
  modl=fltarr(n_chan)
  for i=0,nord do modl=modl+coef(i)*alog(ener)^i

;more documentation type stuff

  if(pflag)then begin
    if(!d.name EQ 'X')then window,2,xsize=512,ysize=512
    if(!d.name EQ 'PS')then dops,'p',1,10.9,10.1
    !p.multi=0
    !mtitle='!17'
    !x.title='!17Energy (keV)'
    !y.title='!17Rate (counts/s/channel)'

    plot_oo,ener(plc_al),objc(plc_al),psym=10
    oplot,ener(plc_al),spec(plc_al),psym=10,color=120
    for i=0,n_elements(plc_al)-1 do oplot,ener(plc_al(i))*[1,1],spec(plc_al(i))+uncr(plc_al(i))*[-1,1],color=120
    oplot,ener(plc_alsi),exp(modl(plc_alsi)),color=220
  for i=0,nels-1 do begin
    oplot,ener(mine(i):maxe(i)),spec(mine(i):maxe(i)),color=60
  endfor

    if(!d.name EQ 'PS')then begin
      device,/close
      spawn,'mv idl.ps back.ps'
    endif
  endif
;finish the bridge

  fillin=fltarr(n_chan)
  fillin(plc_alsi)=exp(modl(plc_alsi))

  spec=spec+fillin

;write to pha file

  ;primary extension
  
  head0=headfits(obj_file,exten=0)
  sxdelpar,head0,['XPROC0','CONTINUE','XDAL0','HISTORY']
  sxaddpar,head0,'OBS_ID','0000000000'
  sxaddpar,head0,'EXP_ID','0000000000000'
  sxaddpar,head0,'REVOLUT',0
  sxaddpar,head0,'OBJECT','Background'
  sxaddpar,head0,'OBSERVER','L.M.A.F.Principe-Kuntz'
  sxaddpar,head0,'RA_OBJ',0.0
  sxaddpar,head0,'DEC_OBJ',0.0
  sxaddpar,head0,'RA_NOM',0.0
  sxaddpar,head0,'DEC_NOM',0.0
  sxaddpar,head0,'REFXCRVL',0.0
  sxaddpar,head0,'REFYCRVL',0.0
  sxaddpar,head0,'RA_PNT',0.0
  sxaddpar,head0,'DEC_PNT',0.0
  sxaddpar,head0,'PA_PNT',0.0
  fxwrite,bck_file,head0

  ;spectrum extension

  chan=indgen(n_elements(spec))
  head1=headfits(obj_file,exten=1)
  ;finish the scaling
  expo=sxpar(head1,'EXPOSURE')
  back=sxpar(head1,'BACKSCAL')
  spec=spec*expo*back
  uncr=uncr*expo*back

  fxbhmake,head2,sxpar(head1,'NAXIS2')
  fxbaddcol,1,head2,chan(0),'CHANNEL'
  fxbaddcol,2,head2,spec(0),'COUNTS'
  fxbaddcol,3,head2,uncr(0),'STAT_ERR'
  
  sxaddpar,head1,'NAXIS1',sxpar(head2,'NAXIS1')
  sxaddpar,head1,'TFIELDS',sxpar(head2,'TFIELDS')
  sxaddpar,head1,'TTYPE1',sxpar(head2,'TTYPE1')
  sxaddpar,head1,'TFORM1',sxpar(head2,'TFORM1')
  sxaddpar,head1,'TTYPE2',sxpar(head2,'TTYPE2')
  sxaddpar,head1,'TFORM2',sxpar(head2,'TFORM2')
  sxaddpar,head1,'TTYPE3',sxpar(head2,'TTYPE3')
  sxaddpar,head1,'TFORM3',sxpar(head2,'TFORM3')
 
  sxdelpar,head1,'SLCTEXPR'
  sxaddpar,head1,'HDUCLAS3','COUNTS'
  sxaddpar,head1,'POISSERR','F'
  sxaddpar,head1,'TUNIT1','channel'
  sxaddpar,head1,'TUNIT2','count'
  sxaddpar,head1,'TUNIT3','count'

  fxbcreate,unit,bck_file,head1
  fxbwritm,unit,['CHANNEL','COUNTS','STAT_ERR'],chan,spec,uncr
  fxbfinish,unit
 
; print,'final stop' 
; stop
  end
