function writeVTK(file_name_output,scalar_values,scalar_label)

#Se escribe los resultados en un archivo vtk

file=open(file_name_output, "w")

write(file, "# vtk DataFile Version 3.0 \n")  #esta linea es la version del archivo y el identificador
write(file, "Unstructured Grid \n")  #esta linea es la version del archivo y el identificador
write(file, "ASCII \n")  #esta linea es la version del archivo y el identificador
write(file, "DATASET UNSTRUCTURED_GRID\n")  #esta linea indica que tipo de estrucutra se incluira. Mallas no estrucutradas para este caso

write(file, "POINTS $Nnodos float \n")  #esta linea inicia indicando los puntos, numero de nodos y el tipo (float, integre)

#A continuaciòn se escriben las coordendas nodales 
for i = 1:Nnodos
    Nodx=  NodalMesh[i,2];
    Nody=  NodalMesh[i,3];
    Nodz=  NodalMesh[i,4];
    write(file, "$Nodx $Nody $Nodz  \n");
end

#Ahora se escribe el conjunto de dataset, en este caso mallas no estructuradas
# su estrucutra es CELLS n size (n= cantidad de elementos, size=tamaño de vertices requeridos)
#Por cada elemento se requieren tres vertices.
S=Nelem*4;
write(file, "CELLS $Nelem $S \n")  #esta linea inicia indicando los puntos, numero de nodos y el tipo (float, integre)

#Ahora se escribe la lista de las celdas con el siguiente orden NumeroPuntos, i,j,k,l,...

for i=1:Nelem
   l= ConeMat[i,2]-1;
   m= ConeMat[i,3]-1;
   n= ConeMat[i,4]-1;
   write(file, "3 $l $m $n \n")  #Siempre va 3 por que son elementos triangulares
end

##Seguido se debe escribir el tipo de celda a la que le corresponde cada elemento. EN este caso es 5 para cada celda
write(file, "CELL_TYPES $Nelem \n")

for i=1:Nelem
    write(file, "5  \n")  #Siempre va 3 por que son elementos triangulares
end

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