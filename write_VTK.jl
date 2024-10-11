function writeVTK(file_name_output,Nnodos,NodalMesh,Nelem,ConeMat,TypeElem,scalar_values,scalar_label)
    #Se escribe los resultados en un archivo vtk
    file=open(file_name_output, "w")
    #Se escribe el encabezado
    write(file, "# vtk DataFile Version 3.0 \n")  #esta linea es la version del archivo y el identificador
    write(file, "Unstructured Grid \n")  #esta linea es la version del archivo y el identificador
    write(file, "ASCII \n")  #esta linea es la version del archivo y el identificador
    write(file, "DATASET UNSTRUCTURED_GRID\n")  #esta linea indica que tipo de estrucutra se incluira. Mallas no estrucutradas para este caso

    write(file, "POINTS $Nnodos float \n")  #esta linea inicia indicando los puntos, numero de nodos y el tipo (float, integre)
    #Se escriben las coordendas nodales 
    for i = 1:Nnodos
        Nodx=  NodalMesh[i,2];
        Nody=  NodalMesh[i,3];
        Nodz=  NodalMesh[i,4];
        write(file, "$Nodx $Nody $Nodz  \n");
    end
    #Ahora se escribe el conjunto de dataset, en este caso mallas no estructuradas
    # su estrucutra es CELLS n size (n= cantidad de elementos, size=tamaño de vertices requeridos)
    VTK_type=0
    n_nod=0
    if TypeElem==2  #Triangulo de 3 nodos
        VTK_type=5
        n_nod=3
    elseif TypeElem==9  #Triangulo de 6 nodos
        VTK_type=22
        n_nod=6
    elseif TypeElem==3  #Cuadrilátero de 4 nodos
        VTK_type=9
        n_nod=4
    elseif TypeElem==10  #Cuadrilátero de 9 nodos
        VTK_type=28
        n_nod=9
    elseif TypeElem==16  #Cuadrilátero de 8 nodos
        VTK_type=23
        n_nod=8
    #Se adicionan los demás elementos de alto orden     
    end
    S=Nelem*(n_nod+1);
    write(file, "CELLS $Nelem $S \n")  #esta linea inicia indicando los puntos, numero de nodos y el tipo (float, integre)
    #Ahora se escribe la lista de las celdas con el siguiente orden NumeroPuntos, i,j,k,l,...
    for i=1:Nelem
    #Se escribe el número de nodos del elemento
    write(file, "$n_nod ")
    for j in 1:n_nod-1
        texto=ConeMat[i,j+1]-1;
        write(file, "$texto ")
    end
    texto=ConeMat[i,n_nod+1]-1;
    write(file, "$texto \n")  #Siempre va 3 por que son elementos triangulares
    end

    #Se escribe el tipo de celda VTK a la que le corresponde cada elemento.
    write(file, "CELL_TYPES $Nelem \n")
    for i=1:Nelem
        write(file, "$VTK_type  \n") 
    end
    #Se escriben los campos escalares asociados a cada nodo
    write(file, "POINT_DATA $Nnodos \n")
    n_scalar_values=size(scalar_values,2)
    for i in 1:n_scalar_values
        labels=scalar_label[i]
        write(file, "SCALARS $labels float \n")
        write(file, "LOOKUP_TABLE default \n")
        for j = 1:Nnodos
            val=  scalar_values[j,i]
            write(file, "$val \n");
        end
    end
    close(file);
end