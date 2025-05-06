using TensorKit

include("../groups.jl")
include("groupelements.jl")


struct H³{G<:Group}
    ω::Function
end

# struct VecGω{G, ω} where {G<:Group, ω<:ThreeCocycle{G}}
#     g::G
# end


