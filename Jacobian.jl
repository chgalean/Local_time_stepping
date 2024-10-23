function Jacobian(chi,eta,e,ConeMat,NodalMesh)
    x,y = nodal_coord(e,ConeMat,NodalMesh)
    n_nod_p=size(x,2)
    #Se definen las derivadas de las funciones base
    dN_dchi=zeros(1,n_nod_p)
    dN_deta=zeros(1,n_nod_p)
    #Se plantea una formulación isoparamétrica para presión y subparamétrica para la velocidad,
    #por lo que el Jacobiano se calcula sobre la base de un elemento QUAD4
    for i in 1:n_nod_p
        _,dN_dchi[1,i],dN_deta[1,i]=N_dN_p(chi,eta,i)
    end
    #Se calculan los términos de la matriz Jacobiana   
    dx_dchi=0.0
    dx_deta =0.0
    dy_dchi=0.0
    dy_deta =0.0
    for i in 1:n_nod_p
        dx_dchi += x[1,i]*dN_dchi[1,i] 
        dx_deta += x[1,i]*dN_deta[1,i]
        dy_dchi += y[1,i]*dN_dchi[1,i]
        dy_deta += y[1,i]*dN_deta[1,i]
    end
    J=[dx_dchi dy_dchi; dx_deta dy_deta]
    detJ=dx_dchi* dy_deta-dy_dchi*dx_deta
    return J,detJ
end
  