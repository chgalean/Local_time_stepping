function K_diff(NodalMesh, ConeMat, e, coef_k, nq)

n_nod=ConeMat[e,1] #Número de nodos del elemento
Kdiff=zeros(n_nod,n_nod)
for l in 1:n_nod
    for m in 1:n_nod
       Kdiff[l,m]=k_lm(l,m,e,ConeMat,NodalMesh,coef_k,nq)
    end
end

return Kdiff;

end
