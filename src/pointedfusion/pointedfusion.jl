using TensorKit
using Revise
include("../groups.jl")
include("groupelements.jl")

abstract type FusionCategory end
abstract type VecGω{G,ω} <: FusionCategory end

abstract type CohomologyGroup{N,G} <: AbelianGroup end
getGroup(::Type{VecGω{G,ω}}) where {G<:Group,ω} = G
struct Irr{𝒞<:FusionCategory} <: Sector
    value::Any
    function Irr{𝒞}(value) where {𝒞<:FusionCategory}
        if 𝒞 <: VecGω
            G = getGroup(𝒞)
            if value isa GroupElement{G}
                new(value)
            else
                throw(ArgumentError("Irr value must be a GroupElement of type $G"))
            end
        else
            throw(ArgumentError("$𝒞 has not been implemented"))
        end
    end
end
ω(a, b, c) = 1
Irr{VecGω{ℤ{3}, ω}}(GroupElement{ℤ{3}}(1))


