function F_l(l,e,ConeMat,NodalMesh,nq)
    #Se definen los puntos y pesos de la cuadratura de Gauss 
    chi_gauss,eta_gauss,pesos=Gauss_qpoints(nq)
    #Se toman las coordenadas nodales  
    x,y = nodal_coord(e,ConeMat,NodalMesh)
    #Se determina el número de nodos en ese elemento
    n_nod=ConeMat[e,1]
    fl=0.0
    for i in 1:nq
        x_coord=0.0
        y_coord=0.0
        for j in 1:n_nod
           N,_,_=N_dN(chi_gauss[i],eta_gauss[i],j)
           x_coord += x[j]*N
           y_coord += y[j]*N
        end
        Q=source_fcn(x_coord,y_coord)
        Nl,_,_=N_dN(chi_gauss[i],eta_gauss[i],l)
        _,detJ=Jacobian(chi_gauss[i],eta_gauss[i],e,ConeMat,NodalMesh)
        fl += pesos[i]*Nl*Q*detJ
    end
    fl *= 0.5
    return fl
end