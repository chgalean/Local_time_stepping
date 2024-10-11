#Este codigo soluciona la ecuaciòn de difusión-advección

#       ∇̇.(-k∇ϕ)+v.∇ϕ= Q     

# por el mètodo de elementos fìnitos utilizando un espacio de elementos triangulares.
# Autor: Cristian Felipe Morales Suàrez
#########################################################################################
#ESPACIO PARA EL LLAMADO DE FUNCIONES Y PAQUETES REQUERIDOS PARA LA SOLUCIÒN DEL SISTEMA
using Plots
using DelimitedFiles
using SparseArrays, LinearAlgebra
include("mesh_import_MSH2.jl")  #Funcion para importar la malla en formato MSH2
include("nodal_coord.jl")       #Funciòn para determinar las coordenadas nodales de un elemento 
include("N_dN.jl")              #Funciòn para calcular las funciones base y sus derivadas 
include("Jacobian.jl")          #Funciòn para calcular el Jacobiano 
include("grad_N.jl")            #Funciòn para calcular el gradiente de una función base 
include("Gauss_qpoints.jl")     #Funciòn para definir los puntos y pesos de la cuadratura de Gauss
include("klm_diff.jl")          #Funciòn para calcuar el componente difusivo de la matriz de rigidez elemental
include("klm_adv.jl")           #Funciòn para calcuar el componente advectivo de la matriz de rigidez elemental
include("K.jl")                 #Funciòn para evaluar la matriz de rigidez global
include("F_l.jl")               #Funciòn para evaluar el vector de cargas elemental
include("F.jl")                 #Funciòn para evaluar el vector de cargas global
include("write_VTK.jl")         #Funciòn para escribir archivos de salida en formato VTK 
include("diff_fcn.jl")          #Funciòn que define el coeficiente de difusión k 
include("source_fcn.jl")        #Funciòn que define el término fuente Q 
include("velocity_fcn.jl")      #Funciòn que define el campo de velocidad advectivo 
#########################################################################################
#PARAMETROS RELACIONADOS AL MODELO
plotmesh_flag=0;  #1 para graficar la malla generada
file_name="Plate_QUAD4_fine"
file_name_mesh=file_name*".msh"
file_name_output=file_name*".vtk"

nq=3;            #Número de puntos de cuadratura a usar en la integración numérica
BC=[0 1;0 0; 0 0; 0 0]  #Se define una matriz con las condiciones de contorno del problema. Cada fila
                        #se refiere a una de los bordes físicos del problema. El valor en la primera columna
                        #define el tipo de condición de borde: 0:Dirichlet 1:Neumann

#######################################################################################
#DISCRETIZACION ESPACIAL
#Se lee el archivo en formato MSH2 que contiene la malla
mesh_file=open(file_name_mesh);
Nnodos,NodalMesh,Nelem,ConeMat,Nfaces,BounCond,TypeElem = mesh_import_MSH2(mesh_file, plotmesh_flag);
#Se crea una matriz de rigidez global y el vector de cargas global
Kglo=spzeros(Nnodos, Nnodos);  #La matriz de rigidez se inicializa como una matriz tipo sparse
Fglo=zeros(Nnodos, 1);

for i in 1:Nelem
    Kele= K(NodalMesh,ConeMat,i,nq)
    Fele=F(NodalMesh,ConeMat,i,nq);
    #Se definen los grados de libertad asociados al elemento
    dofs=[ConeMat[i,2] ConeMat[i,3] ConeMat[i,4] ConeMat[i,5]]
    n_dofs=size(dofs,2)

    #Se realiza el aporte elemental a las matrices globales
    for j in 1:n_dofs
        for k in 1:n_dofs
            Kglo[dofs[j],dofs[k]] += Kele[j,k];
        end
        Fglo[dofs[j]]+= Fele[j];
    end
end
# Se aplican las condiciones de frontera  
# En la matriz de condiciones de frontera la etiqueta 1 de la segunda columna indica que es dirichlet
#Constante de pènalizaciòn
kappa=1e8;
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

#Se escribe el archivo de salida
writeVTK(file_name_output,Nnodos,NodalMesh,Nelem,ConeMat,TypeElem,T,["fi"])
