function N_dN_p(chi,eta,m)
    @assert m <= 4 #Reporta un error si m>4  
    if m==1
      N=0.25*(1-chi)*(1-eta);
      dN_dchi=0.25*(eta-1);
      dN_deta=0.25*(chi-1);
    elseif m==2
      N=0.25*(1+chi)*(1-eta);
      dN_dchi=0.25*(1-eta);
      dN_deta=-0.25*(1+chi);
    elseif m==3
      N=0.25*(1+chi)*(1+eta);
      dN_dchi=0.25*(eta+1);
      dN_deta=0.25*(chi+1);
    elseif m==4
      N=0.25*(1-chi)*(1+eta);
      dN_dchi=-0.25*(1+eta);
      dN_deta=0.25*(1-chi);
    end
    return N,dN_dchi,dN_deta
end