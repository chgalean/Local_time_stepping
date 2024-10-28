function compute_norm(x)
    dim=size(x,2)  #Número de columnas de la matriz x
    long=size(x,1) #Número de filas de la matriz x
    norm_x=zeros(long,1)
    for i in 1:long
      sum=0.0
      for j in 1:dim 
        sum += x[i,j]^2
      end  
      norm_x[i] = sqrt(sum)
    end
    return norm_x
  end