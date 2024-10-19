function N_dN_v(chi,eta,m)
  if m==1
    N=-(1/4)*(1-chi)*(1-eta)*(1+chi+eta)
    dN_dchi=-((eta-1)*(eta+2*chi))/4
    dN_deta=-((chi-1)*(2*eta+chi))/4
  elseif m==2
    N=-(1/4)*(1+chi)*(1-eta)*(1-chi+eta)
    dN_dchi=((eta-1)*(eta-2*chi))/4
    dN_deta=((chi+1)*(2*eta-chi))/4
  elseif m==3
    N=-(1/4)*(1+chi)*(1+eta)*(1-chi-eta)
    dN_dchi=((eta+1)*(eta+2*chi))/4
    dN_deta=((chi+1)*(2*eta+chi))/4
  elseif m==4
    N=-(1/4)*(1-chi)*(1+eta)*(1+chi-eta)
    dN_dchi=-(((eta+1)*(eta-2*chi))/4)
    dN_deta=-(((chi-1)*(2*eta-chi))/4)
  elseif m==5
    N=1/2*(chi^2-1)*(eta-1)
    dN_dchi=chi*(eta-1)
    dN_deta=((chi-1)*(chi+1))/2
  elseif m==6
    N=-(1/2)*(chi+1)*(eta^2-1)
    dN_dchi=-(((eta-1)*(eta+1))/2)
    dN_deta=-((chi+1)*eta)
  elseif m==7
    N=-(1/2)*(chi^2-1)*(eta+1)
    dN_dchi=-(chi*(eta+1))
    dN_deta=-(((chi-1)*(chi+1))/2)
  elseif m==8
    N=1/2*(chi-1)*(eta^2-1)
    dN_dchi=((eta-1)*(eta+1))/2
    dN_deta=(chi-1)*eta
  else
    @assert m < 8
  end
  return N,dN_dchi,dN_deta
end