module Mammo

using Comodo
using FEBio
using Statistics

# Include functions
include("functions.jl")

# Export imported modules for later possible use
export Comodo
export FEBio

# Export functions
export mammodir, breast_surface, tet_breast, hex_breast

end # module Mammo
