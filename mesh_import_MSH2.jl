function mesh_import_MSH2(mesh_file, plotmesh_flag)
    #Se abre el archivo que contiene la informacion de la malla
    MeshInf= readlines(mesh_file)
    close(mesh_file)

    MeshNodes=[]; #Identificador del inicio y final de las coordendas nodales
    MeshElem=[]; #Identificador del inicio y final de de matriz de conexiones elemental
    for (Linea, Info) in enumerate(MeshInf) #En este caso numerate toma un arreglo y devuelve una pareja (indice, valor) del archivo MeshInf con la pareja (indice, valor) se hace una iteraciòn en el ciclo for donde Linea en es indice del archivo y Info recibe el contenido de la linea
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
            return parse(Float32, value)
        end
    end 
    #Primero, cada fila se convierte en un vector fila dependiendo de si 
    VecNodalMesh=[safe_parse.(String.(split(line))) for line in NodalText]'; #La funciòn split divide el texto por separadores. Luego el for line in NodalTex es para que lo haga por cada linea del archivo NodalText

    #Segundo, cada uno de estos vectores debe ser concatenado 
    NodalMesh=vcat(VecNodalMesh...);

    #Se construye la matriz de conectividad
    #la opcion parse convierte una cadena a un numero entero
    NumElem= parse(Int,MeshInf[MeshElem[1]+1]);  

    #Se construye una variable que contenga la informacion unicamente de los elementos
    ElementText=[]
    #Se agrega la informacion de todos los elementos en el siguiente ciclo for
    for i in MeshElem[1]+2:MeshElem[2]-1
        push!(ElementText,MeshInf[i]);
    end

    #Se crea una matriz que irà almacendo el archivo de texto sin tener en cuenta su numeraciòn inicial
    ElementModi=[];  
    #Este ciclo for lee linea por linea, en cada ciclo la linea toma el valor de A y por medio de findfirst encuentra el primer vacio de la cadena
    #luego se agrega al vector ElementModi con la funciòn push! donde el valor a agregar es el indice encontrado por findfirst +1 hasta el final
    n_elememts=length(ElementText)
    for i in 1:n_elememts
    A=ElementText[i]
    IdEmpty=findfirst(isspace,A)
    push!(ElementModi,A[IdEmpty+1:end])
    end 

    BounCond=[]
    ConeMat=[]
    TypeElem=0
    #Se hace un recorrido por todo el listado ElementModi para extraer la información de los elementos y caras externas
    #Los tipos de elementos en formato MSH se encuentran listados en https://gitlab.onelab.info/gmsh/gmsh/blob/master/src/common/GmshDefines.h
    for line in ElementModi
        TypeElem= parse.(Int,(split(line)))[1]
        if TypeElem==1       #Elemento línea de 2 nodos
            n_nods=2
            temp=parse.(Int,(split(line))) 
            lv=size(temp,1)
            temp=[temp[3]; temp[lv-n_nods+1:lv]]
            temp=insert!(temp,1,n_nods)
            temp=transpose(temp)
            push!(BounCond,temp)
        elseif TypeElem==8   #Elemento línea de 3 nodos
            n_nods=3
            temp=parse.(Int,(split(line))) 
            lv=size(temp,1)
            temp=[temp[3]; temp[lv-n_nods+1:lv]]
            temp=insert!(temp,1,n_nods)
            temp=transpose(temp)
            push!(BounCond,temp)
        elseif TypeElem==26   #Elemento línea de 4 nodos
            n_nods=4
            temp=parse.(Int,(split(line))) 
            lv=size(temp,1)
            temp=[temp[3]; temp[lv-n_nods+1:lv]]
            temp=insert!(temp,1,n_nods)
            temp=transpose(temp)
            push!(BounCond,temp)
        elseif TypeElem==27   #Elemento línea de 5 nodos
            n_nods=5
            temp=parse.(Int,(split(line))) 
            lv=size(temp,1)
            temp=[temp[3]; temp[lv-n_nods+1:lv]]
            temp=insert!(temp,1,n_nods)
            temp=transpose(temp)
            push!(BounCond,temp)        
        elseif TypeElem==28   #Elemento línea de 6 nodos
            n_nods=6
            temp=parse.(Int,(split(line))) 
            lv=size(temp,1)
            temp=[temp[3]; temp[lv-n_nods+1:lv]]
            temp=insert!(temp,1,n_nods)
            temp=transpose(temp)
            push!(BounCond,temp)      
        elseif TypeElem==2    #Elemento triangular de 3 nodos
            n_nods=3
            temp=parse.(Int,(split(line))) 
            lv=size(temp,1)
            temp=insert!(temp[lv-n_nods+1:lv],1,n_nods)
            temp=transpose(temp)
            push!(ConeMat,temp); 
        elseif TypeElem==9    #Elemento triangular de 6 nodos
            n_nods=6
            temp=parse.(Int,(split(line))) 
            lv=size(temp,1)
            temp=insert!(temp[lv-n_nods+1:lv],1,n_nods)
            temp=transpose(temp)
            push!(ConeMat,temp); 
        elseif TypeElem==20   #Elemento triangular de 9 nodos
            n_nods=9
            temp=parse.(Int,(split(line))) 
            lv=size(temp,1)
            temp=insert!(temp[lv-n_nods+1:lv],1,n_nods)
            temp=transpose(temp)
            push!(ConeMat,temp);  
        elseif TypeElem==21   #Elemento triangular de 10 nodos
            n_nods=10
            temp=parse.(Int,(split(line))) 
            lv=size(temp,1)
            temp=insert!(temp[lv-n_nods+1:lv],1,n_nods)
            temp=transpose(temp)
            push!(ConeMat,temp); 
        elseif TypeElem==22   #Elemento triangular de 12 nodos
            n_nods=12
            temp=parse.(Int,(split(line))) 
            lv=size(temp,1)
            temp=insert!(temp[lv-n_nods+1:lv],1,n_nods)
            temp=transpose(temp)
            push!(ConeMat,temp); 
        elseif TypeElem==3    #Elemento Cuadrilátero de 4 nodos
            n_nods=4
            temp=parse.(Int,(split(line))) 
            lv=size(temp,1)
            temp=insert!(temp[lv-n_nods+1:lv],1,n_nods)
            temp=transpose(temp)
            push!(ConeMat,temp); 
        elseif TypeElem==10   #Elemento Cuadrilátero de 9 nodos
            n_nods=9
            temp=parse.(Int,(split(line))) 
            lv=size(temp,1)
            temp=insert!(temp[lv-n_nods+1:lv],1,n_nods)
            temp=transpose(temp)
            push!(ConeMat,temp); 
        elseif TypeElem==16   #Elemento Cuadrilátero de 8 nodos
            n_nods=8
            temp=parse.(Int,(split(line))) 
            lv=size(temp,1)
            temp=insert!(temp[lv-n_nods+1:lv],1,n_nods)
            temp=transpose(temp)
            push!(ConeMat,temp); 
        elseif TypeElem==36   #Elemento Cuadrilátero de 16 nodos
            n_nods=16
            temp=parse.(Int,(split(line))) 
            lv=size(temp,1)
            temp=insert!(temp[lv-n_nods+1:lv],1,n_nods)
            temp=transpose(temp)
            push!(ConeMat,temp);    
        elseif TypeElem==37   #Elemento Cuadrilátero de 25 nodos
            n_nods=25
            temp=parse.(Int,(split(line))) 
            lv=size(temp,1)
            temp=insert!(temp[lv-n_nods+1:lv],1,n_nods)
            temp=transpose(temp)
            push!(ConeMat,temp); 
        elseif TypeElem==38   #Elemento Cuadrilátero de 36 nodos
            n_nods=36
            temp=parse.(Int,(split(line))) 
            lv=size(temp,1)
            temp=insert!(temp[lv-n_nods+1:lv],1,n_nods)
            temp=transpose(temp)
            push!(ConeMat,temp); 
        elseif TypeElem==39   #Elemento Cuadrilátero de 12 nodos (incompleto)
            n_nods=12
            temp=parse.(Int,(split(line))) 
            lv=size(temp,1)
            temp=insert!(temp[lv-n_nods+1:lv],1,n_nods)
            temp=transpose(temp)
            push!(ConeMat,temp); 
        elseif TypeElem==40   #Elemento Cuadrilátero de 16 nodos (incompleto)
            n_nods=16
            temp=parse.(Int,(split(line))) 
            lv=size(temp,1)
            temp=insert!(temp[lv-n_nods+1:lv],1,n_nods)
            temp=transpose(temp)
            push!(ConeMat,temp); 
        elseif TypeElem==41   #Elemento Cuadrilátero de 20 nodos (incompleto) 
            n_nods=20
            temp=parse.(Int,(split(line))) 
            lv=size(temp,1)
            temp=insert!(temp[lv-n_nods+1:lv],1,n_nods)
            temp=transpose(temp)
            push!(ConeMat,temp);                                   
        end
    end 
    #Se convierte el arreglo de vectores en una matriz
    BounCond=vcat(BounCond...)
    #Se convierte el arreglo de vectores en una matriz
    ConeMat=vcat(ConeMat...)
    #Se determina el número de elementos
    Nelem= size(ConeMat,1)
    #Se determina el número de caras
    Nfaces= size(BounCond,1)

    if plotmesh_flag==1
        #SE CONSTRUYEN LA MALLA GENERADA
        P2=plot();
        for i in 1: Nelem
            #Se extrae informacion de las coordenadas de los nodos de la forma (x,y) para un nodo dado.
            P= [ConeMat[i,2];ConeMat[i,3];ConeMat[i,4]];
            for j in 1:3
                Nodx=[NodalMesh[P[j],2]; NodalMesh[P[mod(j,3)+1],2]];
                Nody=[NodalMesh[P[j],3]; NodalMesh[P[mod(j,3)+1],3]];
                P2=plot!(Nodx, Nody, legend=false)
            end 
        end

        if plotmesh_flag == 1
        display(P2)
        end
        #savefig("/home/cfmoraless/Codigos/PoissonEquationPlate/malla.png")
    end
    return Nnodos, NodalMesh, Nelem, ConeMat, Nfaces, BounCond, TypeElem;
end