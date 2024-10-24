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
include("N_dN_v.jl")            #Funciòn para calcular las funciones base para la velocidad
include("N_dN_p.jl")            #Funciòn para calcular las funciones base para la presión
include("Jacobian.jl")          #Funciòn para calcular el Jacobiano 
include("grad_N_v.jl")          #Funciòn para calcular el gradiente de una función base de velocidad
include("grad_N_p.jl")          #Funciòn para calcular el gradiente de una función base de presión
include("Gauss_qpoints.jl")     #Funciòn para definir los puntos y pesos de la cuadratura de Gauss
include("Alm_visc.jl")          #Funciòn para calcuar la matriz A elemental
include("A.jl")                 #Funciòn para evaluar la matriz A del término viscoso
include("Blm_imcomp.jl")        #Funciòn para evaluar la matriz B elemental
include("B.jl")                 #Funciòn para evaluar la matriz B del término viscoso
include("F_l.jl")               #Funciòn para evaluar el vector de cargas elemental
include("F.jl")                 #Funciòn para evaluar el vector de cargas global
include("write_VTK.jl")         #Funciòn para escribir archivos de salida en formato VTK 
include("visc_fcn.jl")          #Funciòn que define el coeficiente de difusión k 
include("body_force_fcn.jl")    #Funciòn que define las fuerzas externas sobre el fluido
#########################################################################################
#PARAMETROS RELACIONADOS AL MODELO
plotmesh_flag=0;  #1 para graficar la malla generada
file_name="Plate_QUAD4_coarse"
file_name_mesh=file_name*".msh"
file_name_output=file_name*".vtk"

nq=3;               #Número de puntos de cuadratura a usar en la integración numérica
BC_v=[0 0 0;0 1 0]  #Se define una matriz con las condiciones de contorno de velocidad del problema. Cada fila
                    #se refiere a una de los bordes físicos del problema. El valor en la primera columna
                    #define el tipo de condición de borde: 0:Dirichlet 1:Neumann, la segunda y tercer columna
                    #definen los valores de la velocidad en x y y, respectivamente.

#######################################################################################
#DISCRETIZACION ESPACIAL
#Se lee el archivo en formato MSH2 que contiene la malla
mesh_file=open(file_name_mesh);
Nnodos,NodalMesh,Nelem,ConeMat,Nfaces,BounCond,TypeElem = mesh_import_MSH2(mesh_file, plotmesh_flag);
#Se crea una matriz de rigidez global y el vector de cargas global
Aglo=spzeros(2*Nnodos, 2*Nnodos);  #La matriz A se inicializa como una matriz tipo sparse
Bglo=spzeros(Nnodos, 2*Nnodos);    #La matriz B se inicializa como una matriz tipo sparse
Fglo=zeros(2*Nnodos, 1);

for i in 1:Nelem
    Aele=A(NodalMesh,ConeMat,i,nq)
    Bele=B(NodalMesh,ConeMat,i,nq)
    Fele=F(NodalMesh,ConeMat,i,nq);
    #Se definen los grados de libertad asociados al elemento
    dofs_v=[2*ConeMat[i,2]-1; 
            2*ConeMat[i,2];  
            2*ConeMat[i,3]-1;  
            2*ConeMat[i,3];  
            2*ConeMat[i,4]-1;  
            2*ConeMat[i,4];  
            2*ConeMat[i,5]-1; 
            2*ConeMat[i,5];  
            2*ConeMat[i,6]-1;  
            2*ConeMat[i,6];  
            2*ConeMat[i,7]-1;  
            2*ConeMat[i,7]; 
            2*ConeMat[i,8]-1;  
            2*ConeMat[i,8];  
            2*ConeMat[i,9]-1;  
            2*ConeMat[i,9]]
    n_dofs_v=size(dofs_v,1)
    #Se realiza el aporte elemental de la matriz viscosa elemental Aelem a la matriz global A
    for j in 1:n_dofs_v
        for k in 1:n_dofs_v
            Aglo[dofs_v[j],dofs_v[k]] += Aele[j,k];
        end
        Fglo[dofs_v[j]]+= Fele[j];
    end
    dofs_p=[ConeMat[i,2]; 
            ConeMat[i,3];  
            ConeMat[i,4];  
            ConeMat[i,5]]  
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
#Constante de pènalizaciòn
kappa=1e8;
#Se hace un recorrido por cada una de las caras externas de la malla 
for i in 1:Nfaces
    #Se define el grupo fisico al que pertenece la cara
    phys_grp=BounCond[i,2];
    #Se identifica el tipo de condición de borde correspondiente a ese borde físico
    BC_type=BC_v[phys_grp,1];
    #Se identifica el valor de la condición de borde
    BC_value=BC_v[phys_grp,2:3];
    #Se definen los nodos asociados a la i-esima cara
    nod1=BounCond[i,3]
    nod2=BounCond[i,4]
    nod3=BounCond[i,5]
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
        x= [NodalMesh[nod1,2],NodalMesh[nod3,2]];
        y= [NodalMesh[nod1,3],NodalMesh[nod3,3]];
        #Se calcula la longitud de la cara
        l=sqrt((x[1]-x[2])^2+(y[1]-y[2])^2);
        #Se agrega al vector de cargas
        for j in 1:n_dofs
           #Fglo[dofs[j]]+= 0.5*BC_value*l; 
        end
    end
end 
# Una vez ensambladas las submatrices se ensabla el sistema general y se resuelve
Kglo=[Aglo transpose(Bglo); Bglo zeros(Nnodos,Nnodos)]
display(spy(Kglo)) 
#Fglo=[Fglo;zeros(Nnodos,1)]
T= Aglo\Fglo;
U=T[1:2:2*Nnodos]
V=T[2:2:2*Nnodos]
#C=Bglo*inv(Aglo)
#Bglo*inv(Aglo)*transpose(Bglo)=BA−1F−G,
#T= lu(Kglo) \ Fglo;  #Usando descomposición LU
#T= qr(Kglo) \ Fglo;  #Usando descomposición QR

#Se escribe el archivo de salida
writeVTK(file_name_output,Nnodos,NodalMesh,Nelem,ConeMat,TypeElem,[U V],["U" "V"],[U V],["Velocity"])
