function B(NodalMesh, ConeMat, e, nq)
    n_nod_v=8 #Número de nodos del elemento para interpolar velocidad
    n_nod_p=4 #Número de nodos del elemento para interpolar presión
    Belem=zeros(n_nod_p,2*n_nod_v)
    for l in 1:n_nod_p
       for m in 1:n_nod_v
          Belem[l,2*m-1:2*m] = Blm_imcomp(l,m,e,ConeMat,NodalMesh,nq)
       end
    end
    return Belem;
end