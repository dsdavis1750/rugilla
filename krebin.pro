PRO krebin,oldind,olddep,newind,newdep

npar=n_params()
if(npar EQ 0)then begin
  print,'krebin,oldind,olddep,newind,newdep'
  print,'  program to rebin one vector into another'
  print,'  by flux conservative linear interpolation'
  print,'  inputs:'
  print,'    oldind:old independent vector'
  print,'    olddep:old dependent vector'
  print,'    newind:new independent vector'
  print,'  outputs:'
  print,'    newdep:rebinned dependent vector'
  retall&end

nnew=n_elements(newind)
nold=n_elements(oldind)
newdep=newind*0.0

;check to make sure that both vectors run in the same direction

direc=(oldind(0)-oldind(1))/(newind(0)-newind(1))
if(direc LT 0.0) then begin
  print,'krebin:vectors increase in different directions'
  return
endif

;position 'pointers' such that first new bin's left edge is greater than
;the left edge of the old bin, that is, do the first bin

i=0                                             ;counter for old bins
j=0                                             ;counter for new bins
nlt=newind(0)+(newind(0)-newind(1))/2.0         ;new left (left side of new bin)
olt=oldind(0)+(oldind(0)-oldind(1))/2.0         ;old left (left side of old bin)
if(olt LT nlt)then begin                        ;old extends below new
  while ((olt LT nlt) AND (i LT nold)) do begin ;increment old until new
    i=i+1                                       ;  extends below old
    olt=(oldind(i)+oldind(i-1))/2.0
  endwhile
  if(i EQ nold) then begin
    print,'krebin:vectors are disjoint'
    return
  endif
  i=i-1                                         ;then back up by one bin
  if(i EQ 0)then begin                          ;old left side
    olt=oldind(0)-(oldind(1)-oldind(0))/2.0     ;if first bin
  endif else begin
    olt=(oldind(i)+oldind(i-1))/2.0             ;if not first bin (middle)
  endelse
endif else begin                                ;new extends below old
  while (nlt LT olt) AND (j LT nnew) do begin   ;increment new until old
    j=j+1                                       ;  extends below new
    nlt=(newind(j)+newind(j-1))/2.0
  endwhile
  if(j EQ nnew) then begin
    print,'krebin:vectors are disjoint'
    return
  endif
  j=j-1
  if(j EQ 0)then begin                          ;new left side
    nlt=newind(0)-(newind(1)-newind(0))/2.0     ;if first bin
  endif else begin                              
    nlt=(newind(j)+newind(j-1))/2.0             ;if not first bin (middle)
  endelse
endelse

ort=(oldind(i)+oldind(i+1))/2.0                 ;old right any but last bin

;now do the rest of the newbins

for k=j,nnew-1 do begin                         ;for each new bin
  if(k EQ nnew-1) then begin                    ;new right side
    nrt=newind(k)+(newind(k)-newind(k-1))/2.0   ;if last bin
  endif else begin
    nrt=(newind(k)+newind(k+1))/2.0             ;if middle bin
  endelse
  if(k EQ 0) then begin                         ;new left side
    nlt=newind(0)-(newind(1)-newind(0))/2.0     ;if last bin
  endif else begin
    nlt=(newind(k)+newind(k-1))/2.0             ;if middle bin
  endelse
  if(nrt LT ort)then begin                      ;new completely contained by old
    newdep(k)=olddep(i)*(nrt-nlt)/(ort-olt)     ;add region in common
  endif else begin                              ;new extends above old
    newdep(k)=olddep(i)*(ort-nlt)/(ort-olt)     ;add region in common
    if(i EQ nold-1)then begin
      newdep(k)=newdep(k)*(nrt-nlt)/(ort-nlt)
      return
    endif
    i=i+1                                       ;increment old bin
    olt=(oldind(i)+oldind(i-1))/2.0             ;old left side
    if(i LT nold-1)then begin                   ;old right side
      ort=(oldind(i)+oldind(i+1))/2.0           ;if middle bin
    endif else begin
      ort=oldind(i)+(oldind(i)-oldind(i-1))/2.0 ;if last bin    
    endelse
    while(ort LT nrt) do begin                  ;old contained in new
      if(i EQ nold-1)then begin                 ;last old bin
        newdep(k)=newdep(k)+olddep(i)           ;add all of old bin
        newdep(k)=newdep(k)*(nrt-nlt)/(ort-nlt) ;then normalize
        return 
      endif else begin                          ;not the last bin
        newdep(k)=newdep(k)+olddep(i)           ;add all of old bin
        i=i+1                                   ;increment
        olt=(oldind(i)+oldind(i-1))/2.0         ;old left side
        if(i LT nold-1)then begin               ;old right side
          ort=(oldind(i)+oldind(i+1))/2.0       ;if middle bin
        endif else begin
          ort=oldind(i)+(oldind(i)-oldind(i-1))/2.0   ;if last bin
        endelse
      endelse
    endwhile
    newdep(k)=newdep(k)+olddep(i)*(nrt-olt)/(ort-olt)  ;if old extends above new
  endelse
endfor

return
end
