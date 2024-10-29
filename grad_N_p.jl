function grad_N_p(chi,eta,m,J,detJ)
    @assert m <= 4
    _,dN_dchi,dN_deta=N_dN_p(chi,eta,m)
    grad_N=zeros(1,2)
    grad_N[1,1]=J[2,2]*dN_dchi-J[1,2]*dN_deta
    grad_N[1,2]=-J[2,1]*dN_dchi+J[1,1]*dN_deta
    grad_N = grad_N/detJ
    return grad_N
end
