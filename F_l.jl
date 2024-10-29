function F_l(l,e,ConeMat_P,NodalMesh_P,nq)
  #Se definen los puntos y pesos de la cuadratura en función del número de puntos
  chi_gauss,eta_gauss,pesos=Gauss_qpoints(nq)
  #Se toman las coordenadas nodales  
  x,y = nodal_coord(e,ConeMat_P,NodalMesh_P)
  #Número de nodos de presión
  n_nod_p=4
  fl=zeros(2,1)
  for i in 1:nq
    for k in 1:nq
        x_coord=0.0
        y_coord=0.0
        for j in 1:n_nod_p
            N,_,_=N_dN_p(chi_gauss[i],eta_gauss[k],j)
            x_coord += x[j]*N
            y_coord += y[j]*N
        end
        BF=body_force_fcn(x_coord,y_coord)
        Nl,_,_=N_dN_v(chi_gauss[i],eta_gauss[k],l)
        _,detJ=Jacobian(chi_gauss[i],eta_gauss[k],e,ConeMat_P,NodalMesh_P)
        fl += pesos[i]*pesos[k]*Nl*BF*detJ
    end
  end
  return fl
end