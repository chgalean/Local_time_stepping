#Este codigo soluciona la ecuaciòn de difusión-advección

#       ∇̇.(-k∇ϕ)+v.∇ϕ= Q     

# por el mètodo de elementos fìnitos utilizando un espacio de elementos triangulares.
# Autor: Cristian Felipe Morales Suàrez
#########################################################################################
#ESPACIO PARA EL LLAMADO DE FUNCIONES Y PAQUETES REQUERIDOS PARA LA SOLUCIÒN DEL SISTEMA
using Plots
using DelimitedFiles
using SparseArrays, LinearAlgebra, MUMPS #MKL,  MKL_jll#, MKL, MUMPS, Pardiso,  LinearSolve
include("mesh_import_MSH2.jl")  #Función para importar la malla en formato MSH2
include("nodal_coord.jl")       #Función para determinar las coordenadas nodales de un elemento 
include("N_dN_v.jl")            #Función para calcular las funciones base para la velocidad
include("N_dN_p.jl")            #Función para calcular las funciones base para la presión
include("Jacobian.jl")          #Función para calcular el Jacobiano 
include("grad_N_v.jl")          #Función para calcular el gradiente de una función base de velocidad
include("grad_N_p.jl")          #Función para calcular el gradiente de una función base de presión
include("Gauss_qpoints.jl")     #Función para definir los puntos y pesos de la cuadratura de Gauss
include("Alm_visc.jl")          #Función para calcuar la matriz A elemental
include("A.jl")                 #Función para evaluar la matriz A del término viscoso
include("Blm_imcomp.jl")        #Función para evaluar la matriz B elemental
include("B.jl")                 #Función para evaluar la matriz B del término viscoso
include("F_l.jl")               #Función para evaluar el vector de cargas elemental
include("F.jl")                 #Función para evaluar el vector de cargas global
include("write_VTK.jl")         #Función para escribir archivos de salida en formato VTK 
include("visc_fcn.jl")          #Función que define el coeficiente de difusión k 
include("body_force_fcn.jl")    #Función que define las fuerzas externas sobre el fluido
include("compute_norm.jl")      #Función para el cálculo de la norma de un campo vectorial 
#########################################################################################
#PARAMETROS RELACIONADOS AL MODELO
plotmesh_flag=0;  #1 para graficar la malla generada
file_name="Plate_QUAD4_coarse"
file_name_mesh=file_name*".msh"
file_name_output=file_name*".vtk"

nq=4;               #Número de puntos de cuadratura a usar en la integración numérica
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
print("Finaliza el proceso de ensamble \n")
# Una vez ensambladas las submatrices se ensabla el sistema general y se resuelve
#Kglo=[Aglo transpose(Bglo);Bglo zeros(Nnodos,Nnodos)]
#display(spy(Kglo)) 
############################## MUMPS ##########################################
import MPI
MPI.Init()
root = 0
comm = MPI.COMM_WORLD

mumps0 = MUMPS.Mumps{Float64}(mumps_symmetric, default_icntl, default_cntl64)
if MPI.Comm_rank(comm) == root
    MUMPS.associate_matrix!(mumps0, Aglo)
    MUMPS.associate_rhs!(mumps0, Bglo')
end
MUMPS.factorize!(mumps0)
MUMPS.solve!(mumps0)
MPI.Barrier(comm)
if MPI.Comm_rank(comm) == root
    Schu_comp = Bglo*MUMPS.get_solution(mumps0)
end
finalize(mumps0)

mumps1 = MUMPS.Mumps{Float64}(mumps_symmetric, default_icntl, default_cntl64)
if MPI.Comm_rank(comm) == root
    MUMPS.associate_matrix!(mumps1, Aglo)
    MUMPS.associate_rhs!(mumps1, Fglo)
end
MUMPS.factorize!(mumps1)
MUMPS.solve!(mumps1)
if MPI.Comm_rank(comm) == root
    RHS = Bglo*MUMPS.get_solution(mumps1)
end
finalize(mumps1)

mumps2 = MUMPS.Mumps{Float64}(mumps_symmetric, default_icntl, default_cntl64)
if MPI.Comm_rank(comm) == root
    MUMPS.associate_matrix!(mumps2, Schu_comp)
    MUMPS.associate_rhs!(mumps2, RHS)
end
MUMPS.factorize!(mumps2)
MUMPS.solve!(mumps2)
if MPI.Comm_rank(comm) == root
    p = MUMPS.get_solution(mumps2)
end
finalize(mumps2)

mumps3 = MUMPS.Mumps{Float64}(mumps_symmetric, default_icntl, default_cntl64)
if MPI.Comm_rank(comm) == root
    MUMPS.associate_matrix!(mumps3, Aglo)
    RHS1=Fglo-Bglo'*p
    MUMPS.associate_rhs!(mumps3, RHS1)
end
MUMPS.factorize!(mumps3)
MUMPS.solve!(mumps3)
if MPI.Comm_rank(comm) == root
    UV = MUMPS.get_solution(mumps3)
end
finalize(mumps3)
####
#similar(orig_rhs)
#mumps = MUMPS.Mumps{Float64}(mumps_symmetric)
#Kglo=[Aglo Bglo'; Bglo zeros(Nnodos,Nnodos)]
#MUMPS.associate_matrix!(mumps, Kglo)
#vectx=Vector{Int}(collect(2*Nnodos+1:1:3*Nnodos))
#MUMPS.mumps_schur_complement!(mumps,vectx)

#MUMPS.mumps_solve!(T,Aglo,Fglo) 

#associate_matrix!(mumps, Kglo)
#mumps=MUMPS.SparseArrays(Kglo)
#vectx=spzeros(3*Nnodos,1)
#vectx[1:2*Nnodos,1].=1;
#MUMPS.mumps_schur_complement!(mumps,vectx)
#MPI.Finalize()
############################## MKLPardiso ##########################################
#ps = MKLPardisoSolver()
#set_msglvl!(ps, Pardiso.MESSAGE_LEVEL_ON)
#Kglo=[Aglo Bglo'; Bglo zeros(Nnodos,Nnodos)]
#set_nprocs!(ps, 3)
#Aglo=get_matrix(ps, Aglo, :N)
#BgloT=get_matrix(ps, Bglo, :N)
#T=zeros(2*Nnodos,1)
#RHS=zeros(2*Nnodos,1)
#pardiso(ps, T, Aglo, Fglo)
#S=schur_complement(ps,Kglo,Nnodos)
#S=pardisogetschur(ps)
#S=schur_complement(ps,Kglo,vectx)
#S=pardisogetschur(ps)
#T=solve(ps,Aglo,Bglo'[1:2*Nnodos,1])

#RHS=Bglo*(Aglo\Fglo)
#p=solve(ps,Sch,RHS)
#Kglo=Bglo*Sch
#prob = Aglo\BT[1:2*Nnodos,1]
#T= Aglo\Fglo;
#T= Kglo\Fglo;
#Kglo=[Aglo transpose(Bglo);Bglo zeros(Nnodos,Nnodos)]
#Fglo=[Fglo;zeros(Nnodos,1)]
#T=Kglo\Fglo
print("Finaliza la solución del sistema de ecuaciones \n")

#C=Bglo*inv(Aglo)
#Bglo*inv(Aglo)*transpose(Bglo)=BA−1F−G,
#T= lu(Kglo) \ Fglo;  #Usando descomposición LU
#T= qr(Kglo) \ Fglo;  #Usando descomposición QR

#Se escribe el archivo de salida
U=UV[1:2:2*Nnodos]
V=UV[2:2:2*Nnodos]
norm_V=compute_norm([U V])
writeVTK(file_name_output,Nnodos,NodalMesh,Nelem,ConeMat,TypeElem,norm_V,["|V|"],[U V],["Velocity"])
#writeVTK(file_name_output,Nnodos,NodalMesh,Nelem,ConeMat,TypeElem,[p],["pressure"],[],[])