function assembly(Nnodos_V,Nnodos_P,Nelem_V,ConeMat_V,NodalMesh_P,ConeMat_P,Nfaces_V,BounCond_V,BC_V,nq,kappa)
    #Se crea una matriz de rigidez global y el vector de cargas global
    Aglo=spzeros(2*Nnodos_V, 2*Nnodos_V);  #La matriz A se inicializa como una matriz tipo sparse
    Bglo=spzeros(Nnodos_P, 2*Nnodos_V);    #La matriz B se inicializa como una matriz tipo sparse
    Fglo=spzeros(2*Nnodos_V, 1);

    for i in 1:Nelem_V
        Aele=A(NodalMesh_P,ConeMat_P,i,nq)
        Bele=B(NodalMesh_P,ConeMat_P,i,nq)
        Fele=F(NodalMesh_P,ConeMat_P,i,nq);
        #Se definen los grados de libertad asociados al elemento
        dofs_v=[2*ConeMat_V[i,2]-1; 
                2*ConeMat_V[i,2];  
                2*ConeMat_V[i,3]-1;  
                2*ConeMat_V[i,3];  
                2*ConeMat_V[i,4]-1;  
                2*ConeMat_V[i,4];  
                2*ConeMat_V[i,5]-1; 
                2*ConeMat_V[i,5];  
                2*ConeMat_V[i,6]-1;  
                2*ConeMat_V[i,6];  
                2*ConeMat_V[i,7]-1;  
                2*ConeMat_V[i,7]; 
                2*ConeMat_V[i,8]-1;  
                2*ConeMat_V[i,8];  
                2*ConeMat_V[i,9]-1;  
                2*ConeMat_V[i,9]]
        n_dofs_v=size(dofs_v,1)
        #Se realiza el aporte elemental de la matriz viscosa elemental Aelem a la matriz global A
        for j in 1:n_dofs_v
            for k in 1:n_dofs_v
                Aglo[dofs_v[j],dofs_v[k]] += Aele[j,k];
            end
            Fglo[dofs_v[j]]+= Fele[j];
        end
        dofs_p=[ConeMat_P[i,2]; 
                ConeMat_P[i,3];  
                ConeMat_P[i,4];  
                ConeMat_P[i,5]]  
        n_dofs_p=size(dofs_p,1)
        #Se realiza el aporte elemental de la matriz de incompresibilidad elemental Belem a la matriz global B
        for j in 1:n_dofs_p
            for k in 1:n_dofs_v
                Bglo[dofs_p[j],dofs_v[k]] += Bele[j,k];
            end
        end
    end
    # Se aplican las condiciones de frontera  
    # En la matriz de condiciones de frontera la etiqueta 1 de la segunda columna indica que es dirichlet
    
    #Se hace un recorrido por cada una de las caras externas de la malla 
    for i in 1:Nfaces_V
        #Se define el grupo fisico al que pertenece la cara
        phys_grp=BounCond_V[i,2];
        #Se identifica el tipo de condición de borde correspondiente a ese borde físico
        BC_type=BC_V[phys_grp,1];
        #Se identifica el valor de la condición de borde
        BC_value=BC_V[phys_grp,2:3];
        #Se definen los nodos asociados a la i-esima cara
        nod1=BounCond_V[i,3]
        nod2=BounCond_V[i,4]
        nod3=BounCond_V[i,5]
        #Se definen los grados de libertad asociados a la i-esima cara
        dofs=[2*nod1-1 2*nod1 2*nod2-1 2*nod2 2*nod3-1 2*nod3]
        n_dofs=size(dofs,2)
        if BC_type == 0 #Si se trata de una condición de Dirichlet
            #Se penalizan los grados de libertad asociados a la velocidad en x 
            for j in 1:2:n_dofs
            Aglo[dofs[j], dofs[j]]+= kappa;
            Fglo[dofs[j]]+= BC_value[1]*kappa;
            end
            #Se penalizan los grados de libertad asociados a la velocidad en y 
            for j in 2:2:n_dofs
            Aglo[dofs[j], dofs[j]]+= kappa;
            Fglo[dofs[j]]+= BC_value[2]*kappa;
            end
        else 
        #Se ubican las coordenadas de los nodos que forman parte de la cara
            #x= [NodalMesh[nod1,2],NodalMesh[nod3,2]];
            #y= [NodalMesh[nod1,3],NodalMesh[nod3,3]];
            #Se calcula la longitud de la cara
            #l=sqrt((x[1]-x[2])^2+(y[1]-y[2])^2);
            #Se agrega al vector de cargas
            #for j in 1:n_dofs
            #Fglo[dofs[j]]+= 0.5*BC_value*l; 
            #end
        end
    end 
return Aglo, Bglo, Fglo
end