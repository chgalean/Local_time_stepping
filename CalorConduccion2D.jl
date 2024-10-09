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
using DelimitedFiles
using SparseArrays, LinearAlgebra
include("MeshGenerator0rder_1.jl")  #Funcion para importar la mall
include("K_diff.jl")    #Funciòn para evaluar la matriz de rigidez elemental
include("LoadElemental.jl")         #Funciòn para evaluar el vector de cargas elemental
include("N_dN.jl")                  #Funciòn 
include("Jacobian.jl")              #Funciòn 
include("grad_N.jl")                #Funciòn 
include("k_lm.jl")                  #Funciòn 
include("Write_VTK.jl")                  #Funciòn 
#########################################################################################
#PARAMETROS RELACIONADOS AL MODELO

plotmesh_flag=1;  #1 para graficar la malla generada
file_name="Plate"
file_name_mesh=file_name*".msh"
file_name_output=file_name*".vtk"

coef_k=1.0;    #Coeficiente del término difusivo
Q= 0.0;        # generaciòn de calor
BC=[0 20; 1 50; 0 100; 1 0]    #Se define una matriz con las condiciones de contorno del problema. Cada fila
                        #se refiere a una de los bordes físicos del problema. El valor en la primera columna
                        #define el tipo de condición de borde: 0:Dirichlet 1:Neumann

#######################################################################################
#DISCRETIZACION ESPACIAL

#Mesh es la lectura del archivo .txt que contiene la informaciòn de la malla
mesh_file=open(file_name_mesh);

#Esta funciòn entrega como return entrega:
#Numeros de nodos = Nnodos, 
#Matriz nodal=NodalMesh, 
#Numero de elementos=NumElem,
#Matriz de conectividade=ConeMat


Nnodos, NodalMesh, Nelem, ConeMat, BounCond= MeshGen(mesh_file, plotmesh_flag);
#writedlm("output.txt", BounCond)

#Se crea una matriz de rigidez global y el vector de cargas global
Kglo=spzeros(Nnodos, Nnodos);  #La matriz de rigidez se inicializa como una matriz tipo sparse
Fglo=zeros(Nnodos, 1);

for i in 1:Nelem
    Kele= K_diff(NodalMesh, ConeMat, i, coef_k)
    Fele=LoadElemental(NodalMesh, ConeMat, Q, i);
    #Se definen los grados de libertad asociados al elemento
    dofs=[ConeMat[i,2] ConeMat[i,3] ConeMat[i,4]]
    n_dofs=size(dofs,2)

    #Se realiza el aporte elemental a las matrices globales
    for j in 1:n_dofs
        for k in 1:n_dofs
            Kglo[dofs[j],dofs[k]] += Kele[j,k];
        end
        Fglo[dofs[j]]+= Fele[j];
    end
end

## Se aplican las condiciones de frontera de Dirichlet 
# En la matriz de condiciones de frontera la etiqueta 1 de la segunda columna indica que es dirichlet

#Constante de pènalizaciòn
kappa=1e7;
n_faces=size(BounCond,1)
#Se hace un recorrido por cada una de las caras externas de la malla 
for i in 1:n_faces
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
           Kglo[dofs[j], dofs[j]]+= kappa;
           Fglo[dofs[j]]+= BC_value*kappa;
        end
    else 
       #Se ubican las coordenadas de los nodos que forman parte de la cara
        x= [NodalMesh[nod1,2],NodalMesh[nod2,2]];
        y= [NodalMesh[nod1,3],NodalMesh[nod2,3]];
        #Se calcula la longitud de la cara
        l=sqrt((x[1]-x[2])^2+(y[1]-y[2])^2);

        #Se agrega al vector de cargas
        for j in 1:n_dofs
           Fglo[dofs[j]]+= 0.5*BC_value*l; 
        end
    end

end 

# Una vez acoplado el sistema, se procede a resolver. 

T= Kglo\Fglo;
#T= lu(Kglo) \ Fglo;  #Usando descomposición LU
#T= qr(Kglo) \ Fglo;  #Usando descomposición QR

writeVTK(file_name_output,T,["fi"])




