function StiffnessElemental(NodalMesh, ConeMat, El)

################################################################################
#Se llama la funciòn que calcula el àrea de un elemento

# Se inicia el llenado de la matriz de rigidez elemental

# kel=[ kii kij kik
#       kii kij kik
#       kii kij kik]
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

    a1= (Nk[1]-Nj[1])/(2*AreaEle);
    b1= (Nj[2]-Nk[2])/(2*AreaEle);

    a2= (Ni[1]-Nk[1])/(2*AreaEle);
    b2= (Nk[2]-Ni[2])/(2*AreaEle);

    a3= (Nj[1]-Ni[1])/(2*AreaEle);
    b3= (Ni[2]-Nj[2])/(2*AreaEle);


    K11= b1^2 + a1^2;
    K12= b2*b1 + a2*a1;
    K13= b3*b1 + a3*a1;

    K21= K12;
    K22= b2^2+a2^2;
    K23= b2*b3 + a2*a3;

    K31= K13;
    K32= K23;
    K33= b3^2 + a3^2;

    Kele= AreaEle*[K11 K12 K13; K21 K22 K23;K31 K32 K33];

  
    return Kele;

end
