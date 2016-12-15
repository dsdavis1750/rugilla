PRO aug_plot,rate,urate,hard,uhard,revl,good,c_val

  sls_col,0

;called by get_aug_lim and get_aug_exp

  minr=fix(min(revl)/10.)*10.
  maxr=fix(max(revl)/10.)*10.+10.
  nobs=n_elements(rate)

;plot the rate vs date
 
  !x.range=[minr,maxr]
  !y.range=[0.01,0.06]
  !x.title='!17Revolution'
  !y.title='!17Corner Rate'
 
  plot,revl,rate,psym=6,symsize=0.5
  for i=0,nobs-1 do begin
    oplot,revl(i)*[1,1],rate(i)+urate(i)*[-1,1]
  endfor
  oplot,revl(good),rate(good),psym=6,symsize=0.5,color=60
  if(total(c_val) NE 0)then oplot,c_val(2)*[1,1],c_val(0)*[1,1],psym=1,color=220
 
;plot the hard vs date
 
  !x.range=[minr,maxr]
  !y.range=[0.0,8.0]
  !x.title='!17Revolution'
  !y.title='!17(B-A)/(C-D)'
 
  plot,revl,hard,psym=6,symsize=0.5
  for i=0,nobs-1 do begin
    oplot,revl(i)*[1,1],hard(i)+uhard(i)*[-1,1]
  endfor
  oplot,revl(good),hard(good),psym=6,symsize=0.5,color=60
  if(total(c_val) NE 0)then oplot,c_val(2)*[1,1],c_val(1)*[1,1],psym=1,color=220
 
;plot the hard vs rate
 
  !x.range=[0.01,0.06]
  !y.range=[0.0,8.0]
  !x.title='!17Rate'
  !y.title='!17(B-A)/(C-D)'
 
  plot,rate,hard,psym=6,symsize=0.5
  for i=0,nobs-1 do begin
    oplot,rate(i)+urate(i)*[-1,1],hard(i)*[1,1]
    oplot,rate(i)*[1,1],hard(i)+uhard(i)*[-1,1]
  endfor
  oplot,rate(good),hard(good),psym=6,symsize=0.5,color=60
  if(total(c_val) NE 0)then oplot,c_val(0)*[1,1],c_val(1)*[1,1],psym=1,color=220
 
  !x.range=[0,0]
  !y.range=[0,0]

  end
