function N_dN(chi,eta,m)
  if m==1
    N=chi
    dN_dchi=1.0
    dN_deta=0.0
  elseif m==2
    N=eta
    dN_dchi=0.0
    dN_deta=1.0
  elseif m==3
      N=1-chi-eta
      dN_dchi=-1.0
      dN_deta=-1.0
  else
    @assert m < 3
  end
  return N,dN_dchi,dN_deta
end