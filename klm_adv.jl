function klm_adv(l,m,e,ConeMat,NodalMesh,nq)
    #Se definen los puntos y pesos de la cuadratura de Gauss 
    chi_gauss,eta_gauss,pesos=Gauss_qpoints(nq)
    #Se toman las coordenadas nodales  
    x,y = nodal_coord(e,ConeMat,NodalMesh)
    #Se determina el número de nodos en ese elemento
    n_nod=ConeMat[e,1]
    klm=0.0
    for i in 1:nq
      x_coord=0.0
      y_coord=0.0
      for j in 1:n_nod
         N,_,_=N_dN(chi_gauss[i],eta_gauss[i],j)
         x_coord += x[j]*N
         y_coord += y[j]*N
      end
      v=velocity_fcn(x_coord,y_coord)  
      L_Nm=grad_N(chi_gauss[i],eta_gauss[i],m,e,ConeMat,NodalMesh)
      Nl,_,_=N_dN(chi_gauss[i],eta_gauss[i],l)
      _,detJ=Jacobian(chi_gauss[i],eta_gauss[i],e,ConeMat,NodalMesh)
      klm += pesos[i]*Nl*dot(v,L_Nm)*detJ
    end
    klm *= 0.5
    return klm
end