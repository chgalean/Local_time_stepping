function nodal_coord(e,ConeMat,NodalMesh)
 n_nod=ConeMat[e,1]
 #Se definen los nodos del elemento
 nod=zeros(Int,1,n_nod)
 for i in 1:n_nod
   nod[1,i]=ConeMat[e,1+i] 
 end
 #Se definen las coordenadas de cada uno de los nodos del elemento
 x=zeros(1,n_nod)
 y=zeros(1,n_nod)
 for i in 1:n_nod
    x[1,i]=NodalMesh[nod[1,i],2]
    y[1,i]=NodalMesh[nod[1,i],3]
 end
 return x,y
end