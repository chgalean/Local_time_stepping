function LoadElemental(NodalMesh, ConeMat, Q, El)

    # Se inicia el llenado del vector de cargas elemental

    # Fel=[ fi
    #       fj
    #       fk]
    # Esta matriz es simetrica

    #Se extrae la informaciòn sobre sus respectivos nodos
    NodosEl= ConeMat[El,:];

    #Se obtienen las coordendas de los nodos correspondientes
    Ni=[NodalMesh[NodosEl[2],2], NodalMesh[NodosEl[2],3]];
    Nj=[NodalMesh[NodosEl[3],2], NodalMesh[NodosEl[3],3]];
    Nk=[NodalMesh[NodosEl[4],2], NodalMesh[NodosEl[4],3]];
    
    #Se encuentra el àrea del triàngulo
  
    Dis1=sqrt((Ni[1]-Nj[1])^2 + (Ni[2]-Nj[2])^2);  #Distancias ij
    Dis2=sqrt((Nj[1]-Nk[1])^2 + (Nj[2]-Nk[2])^2);  #Distancias jk
    Dis3=sqrt((Nk[1]-Ni[1])^2 + (Nk[2]-Ni[2])^2);  #Distancias ki

    promdis=(Dis1+Dis2+Dis3)/2;
    AreaEle=sqrt(promdis*(promdis-Dis1)*(promdis-Dis2)*(promdis-Dis3));

    Fele= ((Q*AreaEle)/(3))*[1; 1; 1];

    return Fele
    
end



