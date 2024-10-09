function MeshGen(Mesh, plotmesh)

#using Plots

#se abre el archivo que contiene la informacion de la malla
 
MeshInf= readlines(Mesh)
close(Mesh)

MeshNodes=[]; #Identificador del inicio y final de las coordendas nodales
MeshElem=[]; #Identificador del inicio y final de de matriz de conexiones elemental

for (Linea, Info) in enumerate(MeshInf) #En este caso numerate toma un arreglo y devuelve una pareja (indice, valor) del archivo MeshInf con la pareja (indice, valor) se hace ina iteraciòn en el ciclo for donde Linea en es indice del archivo y Info recibe el contenido de la linea

#Linea por linea del archivo MeshInf Se empezara a buscar donde inicia la coordenadas de los nodos, cuantos nodos hay y donde termina

    if Info ===  "\$Nodes"
        push!(MeshNodes, Linea)
    end

    if Info ===  "\$EndNodes"
        push!(MeshNodes, Linea)
    end

    if Info ===  "\$Elements"
        push!(MeshElem, Linea)
    end

    if Info ===  "\$EndElements"
        push!(MeshElem, Linea)
    end

end

#Se construye una variable que contenga la informacion unicamente de las coordenadas de la malla
NodalText=[]

for i in MeshNodes[1]+2: MeshNodes[2]-1
    push!(NodalText,MeshInf[i]);
end

#Se identifica la cantidad de nodos existentes en la malla 
Nnodos= parse(Int,MeshInf[MeshNodes[1]+1]);  #la opcion parse convierte una cadena a un numero entero

#Una vez identificadas las lineas que contienen las coordenadas que procede a construir la matriz nodal

#Se crea una funciòn que indique si es un valor entero (Int) o Decimal (FLoat64)
function safe_parse(value::String)
    try
        return parse(Int, value)
    catch
        return parse(Float64, value)
    end
end 

#Primero, cada fila se convierte en un vector fila dependiendo de si 
VecNodalMesh=[safe_parse.(String.(split(line))) for line in NodalText]'; #La funciòn split divide el texto por separadores. Luego el for line in NodalTex es para que lo haga por cada linea del archivo NodalText

#Segundo, cada uno de estos vectores debe ser concatenado 
NodalMesh=vcat(VecNodalMesh...);

#AHORA SE CONSTRUYE LA MATRIZ DE CONEXIONES ELEMENTAL

#la opcion parse convierte una cadena a un numero entero
NumElem= parse(Int,MeshInf[MeshElem[1]+1]);  

#Se construye una variable que contenga la informacion unicamente de los elementos
ElementText=[]

#Se agrega la informacion de todos los elementos en el siguiente ciclo for
for i in MeshElem[1]+2: MeshElem[2]-1
    push!(ElementText,MeshInf[i]);
end

#Se crea una matriz que irà almacendo el archivo de texto sin tener en cuenta su numeraciòn inicial
ElementModi=[];  

#Este ciclo for lee linea por linea, en cada ciclo la linea toma el valor de A y por medio de findfirst encuentra el primer vacio de la cadena
#luego se agrega al vector ElementModi con la funciòn push! donde el valor a agregar es el indice encontrado por findfirst +1 hasta el final
for i in 1:length(ElementText)
    A=ElementText[i]
    IdEmpty=findfirst(isspace,A)
    push!(ElementModi,A[IdEmpty+1:end])
end 

ElemBounCon=[]; #matriz que almacenarà los elementos de las fronteras
EleMatrix=[]; #matriz que almacenarà la conectividad de los elementos del dominio

for line in ElementModi
TypeElem= line[1]
    if TypeElem == '1'
        push!(ElemBounCon,line);
    elseif TypeElem == '2'
        push!(EleMatrix,line);
    end
end 

#toma la matriz de numeros en formato texto que se encuentra en ElementText y por medio de split toma linea por linea y los separa por decimales, con esto crea una matriz
BounCond=vcat([parse.(Int,(split(line))) for line in ElemBounCon]'...); 

#toma la matriz de numeros en formato texto que se encuentra en ElementText y por medio de split toma linea por linea y los separa por decimales, con esto crea una matriz
ConeMat=vcat([parse.(Int,(split(line))) for line in EleMatrix]'...) ;

BounCond=BounCond[:,[1,3,5,6]];
ConeMat=ConeMat[:,[1,5,6,7]];

#Se enumeran los elementos del dominio ubicandolos en la columna 1
for i in 1: size(ConeMat,1)
    ConeMat[i,1]= i;
end

#Se enumeran los elementos de la frontera ubicandolos en la columna 1
for i in 1: size(BounCond,1)
    BounCond[i,1]= i;
end

#SE CONSTRUYEN LA MALLA GENERADA
P2=plot();

for i in 1: size(ConeMat,1)

    #Se extrae informacion de las coordenadas de los nodos de la forma (x,y) para un nodo dado.
    P= [ConeMat[i,2];ConeMat[i,3];ConeMat[i,4]];

    for j in 1:3
        Nodx=[NodalMesh[P[j],2]; NodalMesh[P[mod(j,3)+1],2]];
        Nody=[NodalMesh[P[j],3]; NodalMesh[P[mod(j,3)+1],3]];
        P2=plot!(Nodx, Nody, legend=false)
    end 
end

if plotmesh == 1
display(P2)
else

end
#savefig("/home/cfmoraless/Codigos/PoissonEquationPlate/malla.png")

Nelem= ConeMat[end,1]
return Nnodos, NodalMesh, Nelem, ConeMat, BounCond;


end
