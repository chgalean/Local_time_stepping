function Jacobian(chi,eta,e,ConeMat,NodalMesh)
 #Se definen los nodos del elemento
  nod1=ConeMat[e,2]  
  nod2=ConeMat[e,3]
  nod3=ConeMat[e,4]
 #Se definen las coordenadas de cada uno de los nodos del elemento
  x1=NodalMesh[nod1,2]
  y1=NodalMesh[nod1,3]

  x2=NodalMesh[nod2,2]
  y2=NodalMesh[nod2,3]

  x3=NodalMesh[nod3,2]
  y3=NodalMesh[nod3,3]
#Se definen las derivadas de las funciones base
  _,dN1_dchi,dN1_deta=N_dN(chi,eta,1)
  _,dN2_dchi,dN2_deta=N_dN(chi,eta,2)
  _,dN3_dchi,dN3_deta=N_dN(chi,eta,3)
#Se calculan los términos de la matriz Jacobiana   
  dx_dchi = x1*dN1_dchi + x2*dN2_dchi + x3*dN3_dchi 
  dx_deta = x1*dN1_deta + x2*dN2_deta + x3*dN3_deta 

  dy_dchi = y1*dN1_dchi + y2*dN2_dchi + y3*dN3_dchi 
  dy_deta = y1*dN1_deta + y2*dN2_deta + y3*dN3_deta 

  J=[dx_dchi dy_dchi; dx_deta dy_deta]
  detJ=dx_dchi* dy_deta-dy_dchi*dx_deta

  return J,detJ
end
  