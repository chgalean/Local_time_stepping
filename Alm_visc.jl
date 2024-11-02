function Alm_visc(l,m,e,ConeMat_P,NodalMesh_P,nq)
  @assert l <= 8 #Reporta un error si l>8  
  @assert m <= 8 #Reporta un error si m>8 
  #Se definen los puntos y pesos de la cuadratura en función del número de puntos
  chi_gauss,eta_gauss,pesos=Gauss_qpoints(nq)
  #Se toman las coordenadas nodales  
  x,y = nodal_coord(e,ConeMat_P,NodalMesh_P)
  #Número de nodos de presión
  n_nod_p=4
  Alm=zeros(2,2)
  for i in 1:nq
    for k in 1:nq
      x_coord=0.0
      y_coord=0.0
      for j in 1:n_nod_p #Para la interpolación geométrica se usan 4 nodos (los de presión)
        N,_,_=N_dN_p(chi_gauss[i],eta_gauss[k],j)
        x_coord += x[j]*N
        y_coord += y[j]*N
      end
      coeff_nu=visc_fcn(x_coord,y_coord)
      coeff_kappa=kappa_fcn(x_coord,y_coord)
      D=[2*coeff_nu 0 0; 0 2*coeff_nu 0; 0 0 coeff_nu]
      J,detJ=Jacobian(chi_gauss[i],eta_gauss[k],e,ConeMat_P,NodalMesh_P)
      gradNm=grad_N_v(chi_gauss[i],eta_gauss[k],m,J,detJ)
      gradNl=grad_N_v(chi_gauss[i],eta_gauss[k],l,J,detJ)  
      L_Nm=[gradNm[1] 0; 0 gradNm[2];gradNm[2] gradNm[1]]
      L_Nl=[gradNl[1] 0; 0 gradNl[2];gradNl[2] gradNl[1]]
      Nl,_,_=N_dN_v(chi_gauss[i],eta_gauss[k],l)
      Nm,_,_=N_dN_v(chi_gauss[i],eta_gauss[k],m)
      Alm += pesos[i]*pesos[k]*(L_Nl'*D*L_Nm + coeff_nu*coeff_kappa*Nl*Nm*Matrix(1.0I, 2, 2))*detJ
    end
  end
  return Alm
end