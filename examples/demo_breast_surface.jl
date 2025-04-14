using Mammo
using Mammo.Comodo
using Mammo.Comodo.GeometryBasics
using Random
using Statistics

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

# --------------------------------------------------------------------------------
# Visualization
strokewidth = 1 
cmap = Makie.Categorical(:Spectral) 

F1s,V1s = separate_vertices(F1,V1)
C1s = simplex2vertexdata(F1s,C1)

fig = Figure(size=(1200,800))
ax1 = Axis3(fig[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Surface labelling")
hp1 = poly!(ax1,GeometryBasics.Mesh(V1s,F1s), strokewidth=1, color=C1s, strokecolor=:black, shading = FastShading, transparency=false, colormap=cmap)
Colorbar(fig[1, 2], hp1)
ax2 = Axis3(fig[1, 3], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Skin labelling")
hp2 = poly!(ax2,GeometryBasics.Mesh(V1,F1), strokewidth=1, color=C_skin, strokecolor=:black, shading = FastShading, transparency=false, colormap=cmap)
fig