using Mammo
using Mammo.Comodo
using Mammo.Comodo.GLMakie
using Mammo.Comodo.GeometryBasics
using Mammo.Comodo.Rotations
using LinearAlgebra
using Random
using Statistics
using Mammo.FEBio
using Printf 

# Example geometry for a sphere that is cut so some edges are boundary edges
n = 3 # Number of refinement steps of the geodesic sphere
r = 40.0 # hemisphere radius
r1 = r/2.5
r2 = r/7.0
w = (r1-r2)/20.0; #Height of alveola
h = r2/1.5
sy = 1.0 # y-direction squeeze
gravityShiftPercentage = 0.5; #Percentage shift due to gravity

F1, V1, C1, C_skin = breast_surface(n,r,r1,r2,w,h,sy,gravityShiftPercentage)
pointSpacing = mean(edgelengths(F1,V1))

# Meshing with tetgen
v_region = mean(V1)
vol1 = 3.0* pointSpacing^3 / (6.0*sqrt(2.0)) # volume of theoretical perfect tetrahedron

stringOpt = "paAqY"
E,V,CE,Fb,Cb = tetgenmesh(F1,V1; facetmarkerlist=C1, V_regions=[v_region],region_vol=vol1, stringOpt)

# Creat compression plates
zLevelPlate = 35.0
xLevelPlate = 10.0
width_plate = 1.25*r
height_plate = 20.0
depth_plate = 1.25*(2.0*r)
rFillet = 5.0 
pointSpacingPlate = pointSpacing/2.0
contactOffset = pointSpacing/10.0
zMinBreast = minimum([v[3] for v in V1])
V_p=Point{3,Float64}[[xLevelPlate,0.0,zLevelPlate+height_plate],[xLevelPlate,0.0,zLevelPlate],[xLevelPlate+ width_plate,0.0,zLevelPlate],[xLevelPlate+ width_plate,0.0,zLevelPlate+height_plate]]
V_pc = filletcurve(V_p; rMax=rFillet,  constrain_method = :max, n=100, close_loop = false, eps_level = 1e-6)
V_pc = evenly_space(V_pc, pointSpacingPlate; close_loop = false, spline_order = 4)

n_extrude_direction = Vec{3, Float64}(0.0,1.0,0.0) # Extrusion direction
F2,V2 = extrudecurve(V_pc; extent=depth_plate, direction=:both, n=n_extrude_direction, close_loop=false,face_type=:quad)
invert_faces!(F2)

F3 = deepcopy(F2)
V3 = deepcopy(V2)
V3_x = [v[1] for v in V3]
V3 .-= Point3{Float64}(xLevelPlate,0.0,zLevelPlate)
Q = RotXYZ(0.0,π,0.0)
V3 = [Point{3, Float64}(Q*v) for v ∈ V3] 
V3 .+= Point3{Float64}(xLevelPlate+width_plate,0.0,zMinBreast-contactOffset)

indices_chest_nodes = unique(reduce(vcat,Fb[Cb .== 2]))

Fb_contact = Fb[Cb.==1]

compressionAngle = 45.0 * (pi/180)
plateDisplacement_XYZ = Point{3,Float64}(0.0,0.0,20.0)
Q_compression = RotXYZ(compressionAngle,0.0,0.0)
V2 = [Q_compression*v for v in V2]
V3 = [Q_compression*v for v in V3]
plateDisplacement_XYZ = Q_compression * plateDisplacement_XYZ

# Join node sets and fix indices
F2 = [f.+length(V) for f in F2] # Offset indices 
append!(V,V2) # Add nodes
F3 = [f.+length(V) for f in F3] # Offset indices 
append!(V,V3) # Add nodes

# Contact parameters
contactPenalty = 150.0;
laugon = 0;
minaug = 1;
maxaug = 10;
fric_coeff=0.1;


# FEA control settings
numTimeSteps=10 # Number of time steps desired
max_refs=50 # Max reforms
max_ups=0 # Set to zero to use full-Newton iterations
opt_iter=15 # Optimum number of iterations
max_retries=5 # Maximum number of retires
dtmin=(1.0/numTimeSteps)/100.0 # Minimum time step size
dtmax=1.0/numTimeSteps  # Maximum time step size
symmetric_stiffness=0 
min_residual=1e-18

material_type = "Ogden"#"Ogden unconstrained" #"neo-Hookean"

if material_type == "neo-Hookean"
    E_youngs = 0.2e-3 # MPa
    ν = 0.49
elseif material_type == "Ogden unconstrained"
    bulkModulusFactor = 50.0
    c1 = 0.2e-3/3.0 # Shear modulus like parameter, MPa
    m1 = 2 # Set non-linearity
    k = bulkModulusFactor*c1 # Bulk-modulus like parameter, MPa
elseif material_type == "Ogden"
    bulkModulusFactor = 50.0
    c1 = 0.2e-3/3.0 # Shear modulus like parameter, MPa
    m1 = 2 # Set non-linearity
    k = bulkModulusFactor*c1 # Bulk-modulus like parameter, MPa
end

# --------------------------------------------------------------------------------
# FEBIO 

# Set FEBio exec path or name
const FEBIO_EXEC = "febio4" # FEBio executable

######
# Define file names
saveDir = joinpath(mammodir(),"assets","temp") # Main directory to save FEBio input and output files
if !isdir(saveDir)
    mkdir(saveDir)      
end

filename_FEB = joinpath(saveDir,"febioInputFile_01.feb")   # The main FEBio input file
filename_xplt = joinpath(saveDir,"febioInputFile_01.xplt") # The XPLT file for viewing results in FEBioStudio
filename_log = joinpath(saveDir,"febioInputFile_01_LOG.txt") # The log file featuring the full FEBio terminal output stream
filename_disp = "febioInputFile_01_DISP.txt" # A log file for results saved in same directory as .feb file  e.g. nodal displacements
filename_stress = "febioInputFile_01_STRESS.txt"

######
# Define febio input file XML
doc,febio_spec_node = feb_doc_initialize()

aen(febio_spec_node,"Module"; type = "solid") # Define Module node: <Module type="solid"/>

control_node = aen(febio_spec_node,"Control") # Define Control node: <Control>
    aen(control_node,"analysis","STATIC")               
    aen(control_node,"time_steps",numTimeSteps)
    aen(control_node,"step_size",1.0/numTimeSteps)
    aen(control_node,"plot_zero_state",1)
    aen(control_node,"plot_range",@sprintf("%.2f, %.2f",0,-1))
    aen(control_node,"plot_level","PLOT_MAJOR_ITRS")
    aen(control_node,"plot_stride",1)
    aen(control_node,"output_level","OUTPUT_MAJOR_ITRS")
    aen(control_node,"adaptor_re_solve",1)

time_stepper_node = aen(control_node,"time_stepper"; type = "default")
    aen(time_stepper_node,"max_retries",max_retries)
    aen(time_stepper_node,"opt_iter",opt_iter)
    aen(time_stepper_node,"dtmin",dtmin)
    aen(time_stepper_node,"dtmax",dtmax)
    aen(time_stepper_node,"aggressiveness",0)
    aen(time_stepper_node,"cutback",5e-1)
    aen(time_stepper_node,"dtforce",0)

solver_node = aen(control_node,"solver"; type = "solid")
    aen(solver_node,"symmetric_stiffness",symmetric_stiffness)
    aen(solver_node,"equation_scheme",1)
    aen(solver_node,"equation_order","default")
    aen(solver_node,"optimize_bw",0)
    aen(solver_node,"lstol",9e-1)
    aen(solver_node,"lsmin",1e-2)
    aen(solver_node,"lsiter",5)
    aen(solver_node,"max_refs",max_refs)
    aen(solver_node,"check_zero_diagonal",0)
    aen(solver_node,"zero_diagonal_tol",0)
    aen(solver_node,"force_partition",0)
    aen(solver_node,"reform_each_time_step",1)
    aen(solver_node,"reform_augment",0)
    aen(solver_node,"diverge_reform",1)
    aen(solver_node,"min_residual",min_residual)
    aen(solver_node,"max_residual",0)
    aen(solver_node,"dtol",1e-3)
    aen(solver_node,"etol",1e-2)
    aen(solver_node,"rtol",0)
    aen(solver_node,"rhoi",0)
    aen(solver_node,"alpha",1)
    aen(solver_node,"beta",2.5e-01)
    aen(solver_node,"gamma",5e-01)
    aen(solver_node,"logSolve",0)
    aen(solver_node,"arc_length",0)
    aen(solver_node,"arc_length_scale",0)
qn_method_node = aen(solver_node,"qn_method"; type = "BFGS")
    aen(qn_method_node,"max_ups",max_ups)
    aen(qn_method_node,"max_buffer_size",0)
    aen(qn_method_node,"cycle_buffer",0)
    aen(qn_method_node,"cmax",0)

Globals_node   = aen(febio_spec_node,"Globals")

Constants_node = aen(Globals_node,"Constants")
    aen(Constants_node,"R",8.3140000e-06)
    aen(Constants_node,"T",298)
    aen(Constants_node,"F",9.6485000e-05)

Material_node = aen(febio_spec_node,"Material")

material_node = aen(Material_node,"material"; id = 1, name="Material1", type=material_type)
if material_type == "neo-Hookean"
    aen(material_node,"E",E_youngs)
    aen(material_node,"v",ν)
elseif material_type == "Ogden unconstrained"
    aen(material_node,"c1",c1)
    aen(material_node,"m1",m1)
    aen(material_node,"c2",c1)
    aen(material_node,"m2",-m1)
    aen(material_node,"cp",k)
elseif material_type == "Ogden"
    aen(material_node,"c1",c1)
    aen(material_node,"m1",m1)
    aen(material_node,"c2",c1)
    aen(material_node,"m2",-m1)
    aen(material_node,"k",k)
end
    
material_node = aen(Material_node,"material"; id = 2, name="Material2", type="rigid body")
    aen(material_node,"density",1.0)
    aen(material_node,"center_of_mass",join([@sprintf("%.16e",x) for x ∈ mean(V2)],','))

material_node = aen(Material_node,"material"; id = 3, name="Material3", type="rigid body")
    aen(material_node,"density",1.0)
    aen(material_node,"center_of_mass",join([@sprintf("%.16e",x) for x ∈ mean(V3)],','))

Mesh_node = aen(febio_spec_node,"Mesh")

    Nodes_node = aen(Mesh_node,"Nodes"; name="nodeSet_all")
    for (i,v) in enumerate(V)        
        aen(Nodes_node,"node", join([@sprintf("%.16e",x) for x ∈ v],','); id = i)     
    end
    
    # Elements
    Elements_node_1 = aen(Mesh_node,"Elements"; name="Part1", type="tet4")
    for (i,e) in enumerate(E)        
        aen(Elements_node_1,"elem",join([@sprintf("%i",j) for j in e], ','); id = i)
    end
    
    Elements_node_2 = aen(Mesh_node,"Elements"; name="Part2", type="quad4")
    for (i,e) in enumerate(F2)        
        aen(Elements_node_2,"elem",join([@sprintf("%i",j) for j in e], ','); id = i+length(E))
    end

    Elements_node_3 = aen(Mesh_node,"Elements"; name="Part3", type="quad4")
    for (i,e) in enumerate(F3)        
        aen(Elements_node_3,"elem",join(map(string, e), ','); id = i)
    end

    # Node sets
    bcSupportList = "bcSupportList"
    aen(Mesh_node,"NodeSet",join([@sprintf("%i",x) for x ∈ indices_chest_nodes],','); name=bcSupportList)

    surfaceName1 = "Skin"
    Surface_node = aen(Mesh_node,"Surface"; name=surfaceName1)
    for (i,e) in enumerate(Fb_contact)        
        aen(Surface_node,"tri3",join([@sprintf("%i",j) for j in e], ','); id = i)
    end

    surfaceName2 = "Plate1"
    Surface_node = aen(Mesh_node,"Surface"; name=surfaceName2)
    for (i,e) in enumerate(F2)        
        aen(Surface_node,"quad4",join([@sprintf("%i",j) for j in e], ','); id = i)
    end

    surfaceName3 = "Plate2"
    Surface_node = aen(Mesh_node,"Surface"; name=surfaceName3)
    for (i,e) in enumerate(F3)        
        aen(Surface_node,"quad4",join([@sprintf("%i",j) for j in e], ','); id = i)
    end

    surfacePairName1 = "Skin_plate1"
    SurfacePair_node = aen(Mesh_node,"SurfacePair"; name=surfacePairName1)
        aen(SurfacePair_node,"primary",surfaceName1)
        aen(SurfacePair_node,"secondary",surfaceName2)

    surfacePairName2 = "Skin_plate2"
    SurfacePair_node = aen(Mesh_node,"SurfacePair"; name=surfacePairName2)
        aen(SurfacePair_node,"primary",surfaceName1)
        aen(SurfacePair_node,"secondary",surfaceName3)

MeshDomains_node = aen(febio_spec_node, "MeshDomains")
    aen(MeshDomains_node,"SolidDomain"; mat = "Material1", name="Part1")
    
    shelldomain_node = aen(MeshDomains_node,"ShellDomain"; mat = "Material2", name="Part2")
    # aen(shelldomain_node,"shell_thickness",0.1)

    aen(MeshDomains_node,"ShellDomain"; mat = "Material3", name="Part3")
    
Boundary_node = aen(febio_spec_node, "Boundary")
    bc_node = aen(Boundary_node,"bc"; name="zero_displacement", node_set=bcSupportList, type="zero displacement")
        aen(bc_node,"x_dof",1)
        aen(bc_node,"y_dof",1)
        aen(bc_node,"z_dof",1)

Rigid_node = aen(febio_spec_node, "Rigid")
    rigid_bc = aen(Rigid_node,"rigid_bc"; name="RigidFixRot_plate1", type="rigid_fixed")
        aen(rigid_bc,"rb",2)
        aen(rigid_bc,"Ru_dof",1)
        aen(rigid_bc,"Rv_dof",1)
        aen(rigid_bc,"Rw_dof",1)

    rigid_bc = aen(Rigid_node,"rigid_bc"; name="RigidPrescribe_plate1_X", type="rigid_displacement")
        aen(rigid_bc,"rb",2)
        aen(rigid_bc,"dof","x")
        aen(rigid_bc,"value",-plateDisplacement_XYZ[1]; lc=1)
        aen(rigid_bc,"relative",0)
    
    rigid_bc = aen(Rigid_node,"rigid_bc"; name="RigidPrescribe_plate1_Y", type="rigid_displacement")
        aen(rigid_bc,"rb",2)
        aen(rigid_bc,"dof","y")
        aen(rigid_bc,"value",-plateDisplacement_XYZ[2]; lc=1)
        aen(rigid_bc,"relative",0)

    rigid_bc = aen(Rigid_node,"rigid_bc"; name="RigidPrescribe_plate1_Z", type="rigid_displacement")
        aen(rigid_bc,"rb",2)
        aen(rigid_bc,"dof","z")
        aen(rigid_bc,"value",-plateDisplacement_XYZ[3]; lc=1)
        aen(rigid_bc,"relative",0)

    rigid_bc = aen(Rigid_node,"rigid_bc"; name="RigidFixRot_plate2", type="rigid_fixed")
        aen(rigid_bc,"rb",3)
        aen(rigid_bc,"Ru_dof",1)
        aen(rigid_bc,"Rv_dof",1)
        aen(rigid_bc,"Rw_dof",1)

    rigid_bc = aen(Rigid_node,"rigid_bc"; name="RigidPrescribe_plate2_X", type="rigid_displacement")
        aen(rigid_bc,"rb",3)
        aen(rigid_bc,"dof","x")
        aen(rigid_bc,"value",plateDisplacement_XYZ[1]; lc=1)
        aen(rigid_bc,"relative",0)
    
    rigid_bc = aen(Rigid_node,"rigid_bc"; name="RigidPrescribe_plate2_Y", type="rigid_displacement")
        aen(rigid_bc,"rb",3)
        aen(rigid_bc,"dof","y")
        aen(rigid_bc,"value",plateDisplacement_XYZ[2]; lc=1)
        aen(rigid_bc,"relative",0)

    rigid_bc = aen(Rigid_node,"rigid_bc"; name="RigidPrescribe_plate2_Z", type="rigid_displacement")
        aen(rigid_bc,"rb",3)
        aen(rigid_bc,"dof","z")
        aen(rigid_bc,"value",plateDisplacement_XYZ[3]; lc=1)
        aen(rigid_bc,"relative",0)

    Contact_node = aen(febio_spec_node,"Contact")
        contact_node = aen(Contact_node,"contact"; type="sliding-elastic", surface_pair=surfacePairName1)
            aen(contact_node,"two_pass",0)
            aen(contact_node,"laugon",laugon)
            aen(contact_node,"tolerance",0.2)
            aen(contact_node,"gaptol",0.0)
            aen(contact_node,"minaug",minaug)
            aen(contact_node,"maxaug",maxaug)
            aen(contact_node,"search_tol",0.01)
            aen(contact_node,"search_radius",0.1*r)
            aen(contact_node,"symmetric_stiffness",0)
            aen(contact_node,"auto_penalty",1)
            aen(contact_node,"penalty",contactPenalty)
            aen(contact_node,"fric_coeff",fric_coeff)
    
    contact_node = aen(Contact_node,"contact"; type="sliding-elastic", surface_pair=surfacePairName2)
            aen(contact_node,"two_pass",0)
            aen(contact_node,"laugon",laugon)
            aen(contact_node,"tolerance",0.2)
            aen(contact_node,"gaptol",0.0)
            aen(contact_node,"minaug",minaug)
            aen(contact_node,"maxaug",maxaug)
            aen(contact_node,"search_tol",0.01)
            aen(contact_node,"search_radius",0.1*r)
            aen(contact_node,"symmetric_stiffness",0)
            aen(contact_node,"auto_penalty",1)
            aen(contact_node,"penalty",contactPenalty)
            aen(contact_node,"fric_coeff",fric_coeff)            

    LoadData_node = aen(febio_spec_node,"LoadData")
        load_controller_node = aen(LoadData_node,"load_controller"; id=1, name="LC_1", type="loadcurve")
            aen(load_controller_node,"interpolate","LINEAR")
            points_node = aen(load_controller_node,"points")
                aen(points_node,"pt",@sprintf("%.2f, %.2f",0,0))
                aen(points_node,"pt",@sprintf("%.2f, %.2f",1,1))
    
    Output_node = aen(febio_spec_node,"Output")
    
    plotfile_node = aen(Output_node,"plotfile"; type="febio")
        aen(plotfile_node,"var"; type="displacement")
        aen(plotfile_node,"var"; type="stress")
        aen(plotfile_node,"var"; type="relative volume")
        aen(plotfile_node,"var"; type="reaction forces")
        aen(plotfile_node,"var"; type="contact pressure")
        aen(plotfile_node,"compression",@sprintf("%i",0))
    
    logfile_node = aen(Output_node,"logfile"; file=filename_log)
        aen(logfile_node,"node_data"; data="ux;uy;uz", delim=",", file=filename_disp)
        aen(logfile_node,"element_data"; data="s1;s2;s3", delim=",", file=filename_stress)
    # <logfile file="tempModel.txt">
    #   <node_data data="ux;uy;uz" delim="," file="tempModel_disp_out.txt">1, 2, 3, 4, 5, 6, 7, 8, 
    

#######
# Write FEB file
XML.write(filename_FEB, doc)

#######
# Run FEBio
run_febio(filename_FEB,FEBIO_EXEC)

#######
# Import results
DD_disp = read_logfile(joinpath(saveDir,filename_disp))
DD_stress = read_logfile(joinpath(saveDir,filename_stress))
numInc = length(DD_disp)
incRange = 0:1:numInc-1

# Create time varying coordinate vector
VT = Vector{Vector{Point{3,Float64}}}()
@inbounds for i in 0:1:numInc-1
    push!(VT,V .+ [Point{3,Float64}(v) for v in DD_disp[i].data])
end

#######
# --------------------------------------------------------------------------------
# Visualization
strokewidth = 1 
cmap = Makie.Categorical(:Spectral) 

F = element2faces(E) # Triangular faces
CE_F = repeat(CE,inner=4)

Fs,Vs = separate_vertices(F,V)
CE_Vs = simplex2vertexdata(Fs,CE_F)
M = GeometryBasics.Mesh(Vs,Fs)

Fbs,Vbs = separate_vertices(Fb,V)
Cbs = simplex2vertexdata(Fbs,Cb)

fig = Figure(size=(1200,1200))
ax1 = Axis3(fig[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "face_type=:tri, n=0")
hp2 = poly!(ax1,GeometryBasics.Mesh(Vbs,Fbs), strokewidth=1,color=Cbs, strokecolor=:black, shading = FastShading, transparency=false, colormap=cmap)
# hp4 = scatter!(ax1, V1[indVertices_C1_1],markersize=25,color=:red)
Colorbar(fig[1, 2], hp2)

ax2 = Axis3(fig[1, 3], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Cut mesh")
hp2 = poly!(ax2,M, color=CE_Vs, shading = FastShading, transparency=false,strokecolor=:black,strokewidth=strokewidth, overdraw=false,colorrange = (1,2),colormap=cmap)
hp3 = scatter!(ax2, V_p,markersize=15,color=:red)
hp4 = scatter!(ax2, V_pc,markersize=10,color=:black)
hp5 = lines!(ax2, V_pc,linewidth=2,color=:black)
hp6 = poly!(ax2,GeometryBasics.Mesh(V,F2), strokewidth=1,color=:blue, strokecolor=:black, shading = FastShading, transparency=false)
hp7 = poly!(ax2,GeometryBasics.Mesh(V,F3), strokewidth=1,color=:green, strokecolor=:black, shading = FastShading, transparency=false)
# normalplot(ax2,F2,V2)
hp8 = scatter!(ax2, V[indices_chest_nodes],markersize=10,color=:black)

VE  = simplexcenter(E,V)
ZE = [v[3] for v in VE]
Z = [v[3] for v in V]
zMax = maximum(Z)
zMin = minimum(Z)
numSlicerSteps = 3*ceil(Int,(zMax-zMin)/mean(edgelengths(F,V)))

stepRange = range(zMin,zMax,numSlicerSteps)
hSlider = Slider(fig[2, 3], range = stepRange, startvalue = mean(stepRange),linewidth=30)

on(hSlider.value) do z 

    B = ZE .<= z
    indShow = findall(B)
    if isempty(indShow)
        hp2.visible=false        
    else        
        hp2.visible=true
        Fs = element2faces(E[indShow])
        Cs = repeat(CE[indShow],inner=4)
        
        indB = boundaryfaceindices(Fs)        
        Fs = Fs[indB]
        Cs = Cs[indB]
        Fs,Vs = separate_vertices(Fs,V)
        CE_Vs = simplex2vertexdata(Fs,Cs)
        Ms = GeometryBasics.Mesh(Vs,Fs)
        hp2[1] = Ms
        hp2.color = CE_Vs
    end

end
# hSlider.selected_index[]+=1
slidercontrol(hSlider,ax2)

hSlider = Slider(fig[2,4], range = incRange, startvalue = numInc-1,linewidth=30)

nodalColor = lift(hSlider.value) do stepIndex
    norm.(DD_disp[stepIndex].data)
end

M = lift(hSlider.value) do stepIndex    
    return GeometryBasics.Mesh(V.+DD_disp[stepIndex].data,Fb)
end

titleString = lift(hSlider.value) do stepIndex
  "Step: "*string(stepIndex)
end

ax3=Axis3(fig[1,4], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = titleString)

min_p = minp([minp(V) for V in VT])
max_p = maxp([maxp(V) for V in VT])

limits!(ax3, (min_p[1],max_p[1]), 
            (min_p[2],max_p[2]), 
            (min_p[3],max_p[3]))

hp=poly!(ax3,M, strokewidth=2,color=nodalColor, transparency=false, shading = FastShading, colormap = Reverse(:Spectral))#,colorrange=(0,sqrt(sum(displacement_prescribed.^2))))
Colorbar(fig[1,5],hp.plots[1],label = "Displacement magnitude [mm]") 

slidercontrol(hSlider,ax3)

fig
