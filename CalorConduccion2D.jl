#Este codigo soluciona la ecuaciòn de calor dada por

#       ∇̇(∇ϕ)= Q     

# por el mètodo de elementos fìnitos utilizando un espacio de elementos triangulares.
# El dominio discretizado corresponde a una placa con geometrìa rectangular
#
#               \phi=100
#           ________________
#          |                |
#          |                |
#    qn=0  |                | \phi=20
#          |                |
#          |                |
#          |________________|
#                 qn=0
#
# Autor: Cristian Felipe Morales Suàrez
#########################################################################################
#ESPACIO PARA EL LLAMADO DE FUNCIONES Y PAQUETES REQUERIDOS PARA LA SOLUCIÒN DEL SISTEMA

using Plots
include("MeshGenerator0rder_1.jl")  #Funcion para importar la mall
include("StiffnessElemental.jl")    #Funciòn para evaluar la matriz de rigidez elemental
include("LoadElemental.jl")         #Funciòn para evaluar el vector de cargas elemental

#########################################################################################
#PARAMETROS RELACIONADOS AL MODELO

plotmesh=1;  #1 para graficar la malla generada
Q= 0;        # generaciòn de calor
BounDir= 200; #condicion de frontera de Dirichlet
qn= 0;       #Flujo de calor sobre el borde de Neumann


#######################################################################################
#DISCRETIZACION ESPACIAL

#Mesh es la lectura del archivo .txt que contiene la informaciòn de la malla
Mesh=open("/home/cfmoraless/Codigos/PoissonEquationPlate/Plate.msh");

#Esta funciòn entrega como return entrega:
#Numeros de nodos = Nnodos, 
#Matriz nodal=NodalMesh, 
#Numero de elementos=NumElem,
#Matriz de conectividade=ConeMat


Nnodos, NodalMesh, Nelem, ConeMat, BounCond= MeshGen(Mesh, plotmesh);

#Se crea una matriz de rigidez global y el vector de cargas global
Kglo=zeros(Nnodos, Nnodos);
Fglo=zeros(Nnodos, 1);

for i in 1:Nelem

    Kele= StiffnessElemental(NodalMesh, ConeMat, i);
    Fele=LoadElemental(NodalMesh, ConeMat, Q, i);

    dof1=ConeMat[i,2];
    dof2=ConeMat[i,3];
    dof3=ConeMat[i,4];

    Kglo[dof1,dof1] += Kele[1,1];
    Kglo[dof1,dof2] += Kele[1,2];
    Kglo[dof1,dof3] += Kele[1,3];
    
    Kglo[dof2,dof1] += Kele[2,1];
    Kglo[dof2,dof2] += Kele[2,2];
    Kglo[dof2,dof3] += Kele[2,3];

    Kglo[dof3,dof1] += Kele[3,1];
    Kglo[dof3,dof2] += Kele[3,2];
    Kglo[dof3,dof3] += Kele[3,3];

    Fglo[dof1]+= Fele[1];
    Fglo[dof2]+= Fele[2];
    Fglo[dof3]+= Fele[3];

end

## Se aplican las condiciones de frontera de Dirichlet 
# En la matriz de condiciones de frontera la etiqueta 1 de la segunda columna indica que es dirichlet

#Por el mètodo de pènalizaciòn
kappa=1e7;

for i= 1 : size(BounCond,1)

    #Se identifica si un elemento de borde el de condiciòn Dirichlet
    Dirichlet=BounCond[i,2];

    #Etiqueta del elemento
    Elemento=i;

    if Dirichlet == 1

        #Se identifican los nodos que forman el elemento de borde
        dof1= BounCond[Elemento, 3];
        dof2= BounCond[Elemento, 4];

        #Se agrega el factor de penalizaciòn en los elementos correspondientes de la matriz de rigidez global
        Kglo[dof1, dof1]+= kappa;
        Kglo[dof2, dof2]+= kappa;

        #Se agrega el factor de penalizaciòn en los elementos correspondientes al vector de cargas
        Fglo[dof1]+= BounDir*kappa;
        Fglo[dof2]+= BounDir*kappa;

    else i==2

        #Se identifican los nodos que forman el elemento de borde
        dof1= BounCond[Elemento, 3]
        dof2= BounCond[Elemento, 4]
        
        #Se ubican las coordenadas de los nodos que forman parte de la frontera de Neumann
        Ni= [NodalMesh[dof1,2],NodalMesh[dof1,3]];
        Nj= [NodalMesh[dof2,2],NodalMesh[dof2,3]];
        
        #Se encuentra la distancia
        l=sqrt((Ni[1]-Nj[1])^2+(Ni[2]-Nj[2])^2);

        #Se agrega al vector de cargas
        Fglo[dof1]+= 0.5*qn*l; 
        Fglo[dof2]+= 0.5*qn*l; 
        
    end

end 

# Una vez acoplado el sistema, se procede a resolver. 

T= Kglo \ Fglo;


#Se escribe los resultados en un archivo vtk

file=open("/home/cfmoraless/Codigos/PoissonEquationPlate/Temperatures.vtk", "w")

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
 write(file, "SCALARS T float \n")
 write(file, "LOOKUP_TABLE default \n")

 for i = 1:Nnodos
    Tem=  T[i]

    write(file, "$Tem \n");
end

close(file);




