function Mlm(l,m,e,ConeMat,NodalMesh,nq)
    #Se definen los puntos y pesos de la cuadratura de Gauss 
    chi_gauss,eta_gauss,pesos=Gauss_qpoints(nq)
    Mlm=0.0
    for i in 1:nq
      Nl,_,_=N_dN(chi_gauss[i],eta_gauss[i],l)
      Nm,_,_=N_dN(chi_gauss[i],eta_gauss[i],m)
      _,detJ=Jacobian(chi_gauss[i],eta_gauss[i],e,ConeMat,NodalMesh)
      Mlm += pesos[i]*Nl*Nm*detJ
    end
    Mlm *= 0.5
    return Mlm
  end