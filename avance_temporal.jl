function avance_temporal(Cglo,Kglo,Fglo,NodalMesh,BC,BounCond,Nfaces,kappa,theta,Told)
    A=Cglo+theta*Kglo
    b=(Cglo-(1-theta)*Kglo)*Told + Fglo
    # Se aplican las condiciones de frontera  
    #Se hace un recorrido por cada una de las caras externas de la malla 
    for i in 1:Nfaces
        #Se define el grupo fisico al que pertenece la cara
        phys_grp=BounCond[i,2];
        #Se identifica el tipo de condición de borde correspondiente a ese borde físico
        BC_type=BC[phys_grp,1];
        #Se identifica el valor de la condición de borde
        BC_value=BC[phys_grp,2];
        #Se definen los nodos asociados a la i-esima cara
        nod1=BounCond[i,3]
        nod2=BounCond[i,4]
        #Se definen los grados de libertad asociados a la i-esima cara
        dofs=[nod1 nod2]
        n_dofs=size(dofs,2)
        if BC_type == 0 #Si se trata de una condición de Dirichlet
            #Se penalizan los grados de libertad asociados a la cara 
            for j in 1:n_dofs
              A[dofs[j], dofs[j]]+= kappa;
              b[dofs[j]]+= BC_value*kappa;
            end
        else 
            #Se ubican las coordenadas de los nodos que forman parte de la cara
            x= [NodalMesh[nod1,2],NodalMesh[nod2,2]];
            y= [NodalMesh[nod1,3],NodalMesh[nod2,3]];
            #Se calcula la longitud de la cara
            l=sqrt((x[1]-x[2])^2+(y[1]-y[2])^2);
            #Se agrega al vector de cargas
            for j in 1:n_dofs
              b[dofs[j]]+= 0.5*BC_value*l; 
            end
        end
    end 
    Tnew=A\b   
    return Tnew
end