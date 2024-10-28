function Blm_imcomp(l,m,e,ConeMat,NodalMesh,nq)
    #Se definen los puntos y pesos de la cuadratura en función del número de puntos
    chi_gauss,eta_gauss,pesos=Gauss_qpoints(nq)
    Blm=zeros(1,2)
    for i in 1:nq
      for k in 1:nq
        gradNm=grad_N_v(chi_gauss[i],eta_gauss[k],m,e,ConeMat,NodalMesh)
        _,detJ=Jacobian(chi_gauss[i],eta_gauss[k],e,ConeMat,NodalMesh)
        Nl,_,_=N_dN_p(chi_gauss[i],eta_gauss[k],l)
        Blm += -pesos[i]*pesos[k]*Nl*gradNm*detJ
      end
    end
    return Blm
  end