function Gauss_qpoints(nq)
    chi_gauss=[]
    eta_gauss=[]
    pesos=[]
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
    return chi_gauss,eta_gauss,pesos
end