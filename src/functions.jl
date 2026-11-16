# Call required packages
using Comodo
using Comodo.GeometryBasics
using Comodo.Rotations
using FEBio
using Statistics

# Functions 
"""
    mammodir()

Returns package folder 

# Description 
This function simply returns the string for the Mammo path. This is helpful for instance to load items, such as meshes, from the `assets`` folder. 
"""
function mammodir()
    joinpath(@__DIR__, "..")
end

"""
    breast_surface(n::Int,r::Float64, r1=r/7.0,r2=r/7.0,w=(r1-r2)/20.0,h=r2/1.5,sy=1.0,gravityShiftPercentage=0.5)

Creates a parametric breast surface

# Description
This function creates the faces and vertices for a parameterised breast surface model. 
"""
function breast_surface(n::Int,r::Float64, r1=r/7.0,r2=r/7.0,w=(r1-r2)/20.0,h=r2/1.5,sy=1.0,gravityShiftPercentage=0.5)
    rm = mean([r1;r2]) #radius between outer alveola radius and outer nipple radius

    dx = r*gravityShiftPercentage; #X-shift due to gravity
    # Initial hemisphere

    F1,V1,C1 = hemisphere(n,r; face_type=:tri,closed=true)
    
    # Calculate deviation for hemisphere shape
    V1 = [Point(p[1], p[2] * sy, p[3]) for p in V1] # Scale in Y-direction
    indicesTopNodes = unique(reduce(vcat,F1[C1.==1]) )
    dt = [sqrt(v[1]^2+v[2]^2) for v in V1[indicesTopNodes]] # Create distance metric (in 2D) to help define alveola and nipple        

    ge = [p[3]/r for p in V1]  # parameterisation [0-1] 
    V1 = [Point(p[1] + dx * ge[i], p[2] , p[3]) for (i, p) in enumerate(V1)]

    # Create skin vertex labelling for alveola and nipple
    C_skin = zeros(length(V1)) # Vertex color label vector 
    for (j,i) in enumerate(indicesTopNodes)
        #dt = sqrt(V1[i][1]^2+V1[i][2]^2) # Create distance metric (in 2D) to help define alveola and nipple            
        if dt[j] <= r1 && dt[j] > r2 # alveola label
            C_skin[i] = 1               
            g = abs((abs(dt[j] - rm)^3/(rm-r2)^3) - 1.0) # parameterisation [0-1]                
            V1[i] = Point(V1[i][1], V1[i][2], V1[i][3] + w*g)
        elseif dt[j] < r2 # Nipple label        
            C_skin[i] = 2                 
            g = abs((dt[j]^3/r2^3) - 1.0) # parameterisation [0-1]                
            V1[i] = Point(V1[i][1], V1[i][2], V1[i][3] + h*g)
        end 
    end

    Q = RotXYZ(0.0,0.5*pi,0.0)
    V1 = [Point{3, Float64}(Q*v) for v ∈ V1] 
    return F1,V1,C1,C_skin
end

function tet_breast(n,r,r1,r2,w,h,sy,gravityShiftPercentage,volFactor; element_type=Tet4{Int}, stringOpt = "paAqY")
    F1, V1, C1, _ = breast_surface(n,r,r1,r2,w,h,sy,gravityShiftPercentage)
    pointSpacing = mean(edgelengths(F1,V1))

    # Meshing with tetgen
    v_region = mean(V1)
    vol1 = volFactor*pointSpacing^3 / (6.0*sqrt(2.0)) # volume of theoretical perfect tetrahedron

    E, V, _, Fb, Cb = tetgenmesh(F1,V1; facetmarkerlist=C1, V_regions=[v_region],region_vol=vol1, stringOpt=stringOpt, element_type=element_type)
    return E, V, Fb, Cb
end

function hex_breast(r, r1, r2, w, h, sy, gravityShiftPercentage, pointSpacing; f=0.75, nSmooth=25, element_type=Hex8{Int})   
    E, V, F, Fb, Cb = hexhemisphere(r, pointSpacing; f=f, nSmooth=nSmooth, element_type=element_type)    
    
    V = [Point(p[1], p[2] * sy, p[3]) for p in V] # Scale in Y-direction
    indicesTopNodes = unique(reduce(vcat, Fb[Cb.==1]) )
    dt = [sqrt(v[1]^2+v[2]^2) for v in V[indicesTopNodes]] # Create distance metric (in 2D) to help define alveola and nipple        

    dx = r*gravityShiftPercentage
    ge = [p[3]/r for p in V]  # parameterisation [0-1] 
    V = [Point(p[1] + dx * ge[i], p[2] , p[3]) for (i, p) in enumerate(V)]
    rm = mean([r1;r2]) #radius between outer alveola radius and outer nipple radius
    # Create skin vertex labelling for alveola and nipple
    C_skin = zeros(length(V)) # Vertex color label vector 
    for (j,i) in enumerate(indicesTopNodes)
        #dt = sqrt(V1[i][1]^2+V1[i][2]^2) # Create distance metric (in 2D) to help define alveola and nipple            
        if dt[j] <= r1 && dt[j] > r2 # alveola label
            C_skin[i] = 1               
            g = abs((abs(dt[j] - rm)^3/(rm-r2)^3) - 1.0) # parameterisation [0-1]                
            V[i] = Point(V[i][1], V[i][2], V[i][3] + w*g)
        elseif dt[j] < r2 # Nipple label        
            C_skin[i] = 2                 
            g = abs((dt[j]^3/r2^3) - 1.0) # parameterisation [0-1]                
            V[i] = Point(V[i][1], V[i][2], V[i][3] + h*g)
        end 
    end

    Q = RotXYZ(0.0, 0.5*pi, 0.0)
    V = [Q*v for v in V]

    indConstrain = elements2indices(Fb)
    V = smoothmesh_laplacian(F, V, nSmooth, 0.5; constrained_points = indConstrain)
    return E, V, Fb, Cb
end
