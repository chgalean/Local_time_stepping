# Este codigo soluciona la ecuación de Brinkman, que es una mezcla de la ecuación 
# de Stokes y la ecuación de Darcy. El problema de Brinkman plantea

#       ∇̇.(-mu ∇u)+ ∇p + (mu*kappa) * u = F  sobre Ω 
#       ∇̇.(u) = 0   sobre Ω
#       u=Uo   sobre ∂Ω

# La solución se plantea empleando un campo dual de elementos finitos: Quad8 para interpolar 
# el campo de velocidad y Quad4 para interpolar el campo de presión. Esto con el fin de 
# satisfacer la condición inf-sup o condición LBB (Ladyzhenskaya-Babuska-Brezzi). 
# Autor: Carlos Galeano - Cristian Morales
# Universidad Nacional de Colombia
#########################################################################################
#ESPACIO PARA EL LLAMADO DE FUNCIONES Y PAQUETES REQUERIDOS PARA LA SOLUCIÒN DEL SISTEMA
using Plots
using DelimitedFiles
using SparseArrays, LinearAlgebra
using LinearSolve, MUMPS, Base.Threads, IterativeSolvers #MKL,  MKL_jll#, MKL, MUMPS, Pardiso,  LinearSolve
using Dates 
include("mesh_import_MSH2.jl")               #Función para importar la malla en formato MSH2
include("nodal_coord.jl")                    #Función para determinar las coordenadas nodales de un elemento 
include("N_dN_v.jl")                         #Función para calcular las funciones base para la velocidad
include("N_dN_p.jl")                         #Función para calcular las funciones base para la presión
include("Jacobian.jl")                       #Función para calcular el Jacobiano 
include("grad_N_v.jl")                       #Función para calcular el gradiente de una función base de velocidad
include("grad_N_p.jl")                       #Función para calcular el gradiente de una función base de presión
include("Gauss_qpoints.jl")                  #Función para definir los puntos y pesos de la cuadratura de Gauss
include("assembly.jl")                       #Función para el ensamble de las matrices y vectores del sistema de ecuaciones
include("Alm_visc.jl")                       #Función para calcuar la matriz A elemental
include("A.jl")                              #Función para evaluar la matriz A del término viscoso
include("Blm_imcomp.jl")                     #Función para evaluar la matriz B elemental
include("B.jl")                              #Función para evaluar la matriz B del término viscoso
include("F_l.jl")                            #Función para evaluar el vector de cargas elemental
include("F.jl")                              #Función para evaluar el vector de cargas global
include("local2global.jl")                   #Función encargada de llevar los aportes de cada hilo a la matriz global
include("write_VTK.jl")                      #Función para escribir archivos de salida en formato VTK 
include("visc_fcn.jl")                       #Función que define el coeficiente de difusión k 
include("kappa_fcn.jl")                      #Función que define el coeficiente kappa del problema de Brinkman
include("body_force_fcn.jl")                 #Función que define las fuerzas externas sobre el fluido
include("direct_solver_linear_system.jl")    #Función para resolver el sistema de ecuaciones usando el método directo con MUMPS
include("iterative_solver_linear_system.jl") #Función para resolver el sistema de ecuaciones usando un métodos iterativos
include("compute_norm.jl")                   #Función para el cálculo de la norma de un campo vectorial 
#########################################################################################
#PARAMETROS RELACIONADOS AL MODELO
plotmesh_flag=0;  #1 para graficar la malla generada
file_name="Plate_coarse"
file_name_mesh_P=file_name*"_P.msh"
file_name_mesh_V=file_name*"_V.msh"
file_name_output_P=file_name*"_P.vtk"
file_name_output_V=file_name*"_V.vtk"

nq=3;                            #Número de puntos de cuadratura a usar en la integración numérica
BC_V=[0 0 0;0 1 0;0 0 0; 0 1 0]  #Se define una matriz con las condiciones de contorno de velocidad del problema. Cada fila
                                 #se refiere a una de los bordes físicos del problema. El valor en la primera columna
                                 #define el tipo de condición de borde: 0:Dirichlet 1:Neumann, la segunda y tercer columna
                                 #definen los valores de la velocidad en x y y, respectivamente.
#Constante de pènalizaciòn
kappa=1e12;
#######################################################################################
#DISCRETIZACION ESPACIAL
#Se lee el archivo en formato MSH2 que contiene la malla QUAD8 para interpolar velocidad
mesh_file_V=open(file_name_mesh_V);
Nnodos_V,NodalMesh_V,Nelem_V,ConeMat_V,Nfaces_V,BounCond_V,TypeElem_V = mesh_import_MSH2(mesh_file_V, plotmesh_flag);
#Se lee el archivo en formato MSH2 que contiene la malla QUAD4 para interpolar presión
mesh_file_P=open(file_name_mesh_P);
Nnodos_P,NodalMesh_P,Nelem_P,ConeMat_P,Nfaces_P,BounCond_P,TypeElem_P = mesh_import_MSH2(mesh_file_P, plotmesh_flag);

times=Dates.format(now(), "HH:MM")
print("Inicia el proceso de ensamble:  " * times * "\n")
print("Este proceso usa ", Threads.nthreads(), " hilos","\n")
Aglo, Bglo, Fglo = assembly(Nnodos_V,Nnodos_P,Nelem_V,ConeMat_V,NodalMesh_P,ConeMat_P,Nfaces_V,BounCond_V,BC_V,nq,kappa)
#display(spy([Aglo Bglo'; Bglo zeros(Nnodos_P,Nnodos_P)],title="Sparsity pattern of KG"))
times=Dates.format(now(), "HH:MM")
print("Finaliza el proceso de ensamble:  " * times * " \n")

#Se soluciona el sistema de ecuaciones
times=Dates.format(now(), "HH:MM")
print("Inicia la solución del sistema de ecuaciones:  " * times * " \n")
import MPI
UV, p = direct_solver_linear_system(Aglo, Bglo, Fglo)
#UV, p = iterative_solver_linear_system(Aglo, Bglo, Fglo)
times=Dates.format(now(), "HH:MM")
print("Finaliza la solución del sistema de ecuaciones:  " * times * " \n")

#Se escriben los archivos de salida
U=UV[1:2:2*Nnodos_V]
V=UV[2:2:2*Nnodos_V]
norm_V=compute_norm([U V])
writeVTK(file_name_output_V,Nnodos_V,NodalMesh_V,Nelem_V,ConeMat_V,TypeElem_V,norm_V,["|V|"],[U V],["Velocity"])
writeVTK(file_name_output_P,Nnodos_P,NodalMesh_P,Nelem_P,ConeMat_P,TypeElem_P,p,["pressure"],[],[])
times=Dates.format(now(), "HH:MM")
print("Cálculo terminado:  " * times * " \n")