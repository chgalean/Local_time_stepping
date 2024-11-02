function local2global(Aglo_vector,Bglo_vector,Fglo_vector,Aglo,Bglo,Fglo,n_threads)
    for i in 1:n_threads
       Aglo += Aglo_vector[i]
       Bglo += Bglo_vector[i] 
       Fglo += Fglo_vector[i]
    end
    return Aglo, Bglo, Fglo
end