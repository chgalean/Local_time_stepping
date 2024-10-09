function k_lm(l,m,e,ConeMat,NodalMesh,coef_k,nq)

if nq==1
  chi_gauss=[1/3]
  eta_gauss=[1/3]
  pesos=[1]
elseif nq==3
    chi_gauss=[1/6 1/6 2/3]
    eta_gauss=[1/6 2/3 1/6]
    pesos=[1/3 1/3 1/3]        
elseif nq==4
    chi_gauss=[1/3 1/5 1/5 3/5]
    eta_gauss=[1/3 1/5 3/5 1/5]
    pesos=[-0.5625 0.52083333333333 0.52083333333333 0.52083333333333]
end

klm=0.0
for i in 1:nq
  L_Nm=grad_N(chi_gauss[i],eta_gauss[i],m,e,ConeMat,NodalMesh)
  L_Nl=grad_N(chi_gauss[i],eta_gauss[i],l,e,ConeMat,NodalMesh)  
  _,detJ=Jacobian(chi_gauss[i],eta_gauss[i],e,ConeMat,NodalMesh)
  klm += pesos[i]*dot(L_Nl,L_Nm)*detJ
end
klm *= 0.5*coef_k
return klm
end