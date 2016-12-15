FUNCTION xmm_eer,frac,r_core,alpha

  rc=r_core
  al=alpha

  r=rc*sqrt(((1.-frac*(1.-((1.+((5.*60./rc)^2))^(1.-al))))^(1./(1.-al)))-1)
  return,r
  end
"
