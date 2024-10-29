function Blm_imcomp(l,m,e,ConeMat_P,NodalMesh_P,nq)
    @assert l <= 4 #Reporta un error si l>4  
    @assert m <= 8 #Reporta un error si m>8  
    #Se definen los puntos y pesos de la cuadratura en función del número de puntos
    chi_gauss,eta_gauss,pesos=Gauss_qpoints(nq)
    Blm=zeros(1,2)
    for i in 1:nq
      for k in 1:nq
        J,detJ=Jacobian(chi_gauss[i],eta_gauss[k],e,ConeMat_P,NodalMesh_P)
        gradNm=grad_N_v(chi_gauss[i],eta_gauss[k],m,J,detJ)
        Nl,_,_=N_dN_p(chi_gauss[i],eta_gauss[k],l)
        Blm -= pesos[i]*pesos[k]*Nl*gradNm*detJ
      end
    end
    return Blm
  end