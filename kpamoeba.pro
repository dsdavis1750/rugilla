;name:amoerr.am.pro
;date:22 January 1991
;purp:amoeba fitting routine
;auth:translated from Numerical Recipes into IDL by K.D.Kuntz
;mods:created outter loop structure to avoid local minima
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
PRO load,w,f,e

  COMMON data,indep,dep,edep,res

  indep=w
  dep=f
  edep=e
  end

;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
PRO subamoeba,func,vertex,value,mask,tolr,maxiter,iter

;inputs
;  tolr=the goodness to which the minimized value is to be found
;  maxiter=maximum number of iterations to be executed
;  vertex=array(ndim,ndim+1) of ndim+1 vertices in free variable space
;  value=array(ndim) of value of minimized function at each vertex
;  mask=mask for changeable parameters
;output
;  iter=number of iterations required (-1 if maxed out without a solution)
;  vertex=set of vertices within tolerance

  COMMON data,cindep,cdep,cedep,cres

;hardwired parameters
  alpha=1.0      ;reflection parameter
  beta=0.5       ;contraction parameter
  gamma=2.0      ;expansion parameter

;setup
  ndim=n_elements(value)-1
  mdim=ndim-1
  numvert=ndim+1
  pbar=fltarr(ndim)
  parms=fltarr(ndim)
  newvert=fltarr(ndim)
  newnewvert=fltarr(ndim)
  iter=0
  nsflag=0
  puflag=n_elements(parms)
  aflag=0

  start:
  lo=1
  for j=0,mdim do begin
    parms(j)=total(vertex(j,*))/(ndim+1)
  endfor

  if(value(0) GT value(1))then begin
    hi=0
    nhi=1
  endif else begin
    hi=1
    nhi=0
  endelse

;find point with the highest chisqr
  for i=0,numvert-1 do begin
    if(value(i) LT value(lo))then lo=i
    if(value(i) GT value(hi))then begin
      nhi=hi
      hi=i
    endif else begin
      if(value(i) GT value(nhi))then begin
        if(i NE hi) then nhi = i
      endif
    endelse
  endfor

;compute fractional range and test for return
  if(iter GE maxiter)then begin
    iter=-1
    return
  endif
  rtol=2.0*ABS(value(hi)-value(lo))/(abs(value(hi))+abs(value(lo)))
  if((rtol LT tolr)AND(iter GT 0))then begin
    if(nsflag NE 0) then begin
      return
    endif
  endif
  puflag=ndim
  for i=0,mdim do begin
    if(mask(i) ne 0) then begin
      if(2.0*abs(vertex(i,hi)-pbar(i))/(abs(vertex(i,hi))+abs(pbar(i))) lt 1.0e-8) then begin
        puflag=puflag-1
      endif
    endif
  endfor
  if(puflag LT 1) then begin
    iter=iter*(-1)
    return
  endif
  nsflag=0
  iter=iter+1

;compute average vector (excluding high point)
  pbar=pbar*0.0
  for i=0,numvert-1 do begin
    if(i NE hi) then pbar=pbar+vertex(*,i)
  endfor
  pbar=pbar/ndim
  newvert=mask*alpha*(pbar-vertex(*,hi))+pbar
  newval=call_function(func,newvert)

  if(newval LE value(lo))then begin
    newnewvert=mask*gamma*(newvert-pbar)+pbar
          newnewval=call_function(func,newnewvert)
    if(newnewval LT value(lo))then begin
      vertex(0,hi)=newnewvert
      value(hi)=newnewval
    endif else begin
      vertex(0,hi)=newvert
      value(hi)=newval
    endelse
  endif else begin
    if(newval GE value(nhi))then begin
      if(newval LT value(hi))then begin
        vertex(0,hi)=newvert
        value(hi)=newval
        aflag=1
      endif
      newnewvert=mask*beta*(vertex(*,hi)-pbar)+pbar
      newnewval=call_function(func,newnewvert)
      if(newnewval LT value(hi)) then begin
        vertex(0,hi)=newnewvert
        value(hi)=newnewval
        aflag=0
      endif else begin
        for i=0,numvert-1 do begin
          if(i NE lo) then begin
            nsflag=1
            newvert=0.5*(vertex(*,i)+vertex(*,lo))
            vertex(0,i)=newvert
            value(i)=call_function(func,newvert)
          endif
        endfor
        aflag=0
      endelse
      aflag=0
    endif else begin
      vertex(0,hi)=newvert
      value(hi)=newval
    endelse
  endelse

  goto,start
  return
  end

;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
PRO kpamoeba,func,iparms,dparms,mask,tolr,maxiter,reslo,reshi,iter,chilo,chihi

;inputs
;  iparms=initial parameters for fitting functions
;  dparms=distance along each axis in phase space to be searched
;  mask=mask for disabbling change to a parameter
;        changeable parameters=1, unchangeable parameters=0
;  tolr=(max(dep-fit)-min(dep-fit))/avg(dep-fit)
;  maxiter=max number of iterations to execute
;outputs
;  reslo=set of best parameters
;  reshi=set of worst allowable parameters
;  iter=number of iterations run (-1 if maxed out without hitting tolr)
;  chilo=chisqr for best fit
;  chihi=chisqr for worst allowable fit

  COMMON data,cindep,cdep,cedep,cres

  iter=0
  fact=-1

  ndim=n_elements(iparms)
  vertex=dblarr(ndim,ndim+1)
  value=dblarr(ndim+1)
  res=dblarr(ndim)
  parms=iparms

  chisqr=call_function(func,parms)

  loop:
  vertex(0,0)=parms
  value(0)=call_function(func,parms)
  fact=fact*(-1)
  for i=1,ndim do begin
    vertex(0,i)=parms
    vertex(i-1,i)=vertex(i-1,i)+dparms(i-1)*mask(i-1)*fact
    temp=vertex(*,i)
    value(i)=call_function(func,temp)
  endfor

  subamoeba,func,vertex,value,mask,tolr,maxiter,subiter
  if(subiter EQ -1)then subiter=maxiter
  iter=iter+subiter
  oldchisqr=chisqr
  rank=sort(value)

  res=vertex(*,rank(0))

  chisqr=value(rank(0))
  rtol=abs(2.0*(oldchisqr-chisqr)/(oldchisqr+chisqr))
  ;condition for looping
  if(rtol GT tolr) then begin
    parms=res
    goto,loop
  endif

  reshi=vertex(*,rank(ndim))
  chihi=call_function(func,reshi)
  reslo=res
  chilo=chisqr

 end


