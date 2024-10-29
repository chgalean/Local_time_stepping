function A(NodalMesh_P, ConeMat_P, e, nq)
 n_nod_v=8 #Número de nodos para interpolar la velocidad
 Aelem=zeros(2*n_nod_v,2*n_nod_v)
 for l in 1:n_nod_v
    for m in 1:n_nod_v
       Aelem[2*l-1:2*l , 2*m-1:2*m] = Alm_visc(l,m,e,ConeMat_P,NodalMesh_P,nq)
    end
 end
 return Aelem;
end
