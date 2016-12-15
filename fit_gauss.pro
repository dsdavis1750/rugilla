FUNCTION fit_gauss,parms

Common gauss_com,gauss_ind,gauss_dep

  parms=abs(parms)
  res=parms(0)*exp(-(gauss_ind-parms(1))^2/2./parms(2)/parms(2))
  chi=total((gauss_dep-res)^2)/n_elements(res)
  return,chi
  end
