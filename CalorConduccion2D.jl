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
include("Clm.jl")               #Funciòn para calcuar la matriz capacitiva elemental
include("K.jl")                 #Funciòn para evaluar la matriz de rigidez global
include("C.jl")                 #Funciòn para evaluar la matriz capacitiva global
include("F_l.jl")               #Funciòn para evaluar el vector de cargas elemental
include("F.jl")                 #Funciòn para evaluar el vector de cargas global
include("write_VTK.jl")         #Funciòn para escribir archivos de salida en formato VTK 
include("diff_fcn.jl")          #Funciòn que define el coeficiente de difusión k 
include("source_fcn.jl")        #Funciòn que define el término fuente Q 
include("velocity_fcn.jl")      #Funciòn que define el campo de velocidad advectivo 
include("avance_temporal.jl")   #Funciòn con el esquema de avance temporal
#########################################################################################
#PARAMETROS RELACIONADOS AL MODELO
plotmesh_flag=0;  #1 para graficar la malla generada
file_name="Plate_TRIANG3_coarse"
file_name_mesh=file_name*".msh"

nq=3;            #Número de puntos de cuadratura a usar en la integración numérica
BC=[0 1;0 0; 0 0; 0 0]  #Se define una matriz con las condiciones de contorno del problema. Cada fila
                        #se refiere a una de los bordes físicos del problema. El valor en la primera columna
                        #define el tipo de condición de borde: 0:Dirichlet 1:Neumann
delta_t=0.01   #Paso de tiempo
Ntime_step=200   #Número de pasos de tiempo
save_step=10     #Intervalo de guardado de datos
theta=1.0      #Parámetro para definir el método: 0   =>Forward Euler
                        #                                  0.5 =>Crank-Nicolson
                        #                                  1   =>Backward Euler
kappa=1e9;     #Constante de pènalizaciòn
#######################################################################################
#DISCRETIZACION ESPACIAL
#Se lee el archivo en formato MSH2 que contiene la malla
mesh_file=open(file_name_mesh);
Nnodos,NodalMesh,Nelem,ConeMat,Nfaces,BounCond,TypeElem = mesh_import_MSH2(mesh_file, plotmesh_flag);
#Se crea el vector con las condiciones iniciales
Tini=0.0*ones(Nnodos,1)
#Se crea el primer archivo de salida
file_name_output=joinpath(pwd()*"/Sol",file_name*"_0.vtk") 
writeVTK(file_name_output,Nnodos,NodalMesh,Nelem,ConeMat,TypeElem,Tini,["T"])
#Se crea una matriz de rigidez global y el vector de cargas global
Kglo=spzeros(Nnodos, Nnodos);  #La matriz de rigidez se inicializa como una matriz tipo sparse
Cglo=spzeros(Nnodos, Nnodos);  #La matriz capacitiva se inicializa como una matriz tipo sparse
Fglo=zeros(Nnodos, 1);
#Se hace un recorrido por los elementos para realizar en ensamble de las matrices y vectores globales
for i in 1:Nelem
    Kele=K(NodalMesh,ConeMat,i,nq)
    Cele=C(NodalMesh,ConeMat,i,nq)
    Fele=F(NodalMesh,ConeMat,i,nq);
    #Se definen los grados de libertad asociados al elemento
    dofs=[ConeMat[i,2] ConeMat[i,3] ConeMat[i,4]]
    n_dofs=size(dofs,2)
    #Se realiza el aporte elemental a las matrices globales
    for j in 1:n_dofs
        for k in 1:n_dofs
            Kglo[dofs[j],dofs[k]] += Kele[j,k];
            Cglo[dofs[j],dofs[k]] += Cele[j,k];
        end
        Fglo[dofs[j]]+= Fele[j];
    end
end
#Proceso de avance temporal
Cglo *=(1/delta_t)
global cont=0
global T=Tini
for t in 1:Ntime_step
    global cont +=1
    global T=avance_temporal(Cglo,Kglo,Fglo,NodalMesh,BC,BounCond,Nfaces,kappa,theta,T)
    if cont==10
      #Se escribe el archivo de salida
      file_name_VTK=joinpath(pwd()*"/Sol",file_name*"_"*string(t)*".vtk") 
      writeVTK(file_name_VTK,Nnodos,NodalMesh,Nelem,ConeMat,TypeElem,T,["T"])
      global cont=0
    end
end
#T= lu(Kglo) \ Fglo;  #Usando descomposición LU
#T= qr(Kglo) \ Fglo;  #Usando descomposición QR