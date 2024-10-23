function nodal_coord(e,ConeMat,NodalMesh)
 #Se trabaja una formulación subparamétrica, de modo que para la interpolación geométrica
 #se usan solo 4 ndoos (los de presión)
 n_nod_p=4
 #Se definen los nodos del elemento
 nod=zeros(Int,1,n_nod_p)
 for i in 1:n_nod_p
  #Se usan solo los 4 primeros nodos (los de presión)
   nod[1,i]=ConeMat[e,1+i] 
 end
 #Se definen las coordenadas de cada uno de los nodos del elemento
 x=zeros(1,n_nod_p)
 y=zeros(1,n_nod_p)
 for i in 1:n_nod_p
    x[1,i]=NodalMesh[nod[1,i],2]
    y[1,i]=NodalMesh[nod[1,i],3]
 end
 return x,y
end