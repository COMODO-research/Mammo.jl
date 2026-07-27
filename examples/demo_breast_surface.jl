using Mammo
using Mammo.Comodo
using Mammo.Comodo.GeometryBasics
using Mammo.Comodo.GLMakie

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

F1s, V1s = separate_vertices(F1, V1)
C1s = simplex2vertexdata(F1s, C1)

fig = Figure(size=(1200,800))
ax1 = AxisGeom(fig[1, 1], title = "Surface labelling")
hp1 = meshplot!(ax1, F1s, V1s, strokewidth=1, color=C1s, strokecolor=:black, colormap=cmap)
Colorbar(fig[1, 2], hp1)
ax2 = AxisGeom(fig[1, 3], title = "Skin labelling")
hp2 = meshplot!(ax2, F1, V1, strokewidth=1, color=C_skin, strokecolor=:black, colormap=cmap)
fig