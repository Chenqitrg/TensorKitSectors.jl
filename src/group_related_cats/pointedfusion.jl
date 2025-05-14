struct VecGωIrr{G<:Group, ω} <: Sector
    g::G
    function VecGωIrr{G, ω}(g) where {G<:Group, ω}
        if g isa G
            new{G, ω}(g)
        else
            throw(ArgumentError("The Irr value $g must be a GroupElement of type $G"))
        end
    end
end

rank(::Type{VecGωIrr{G, ω}}) where {G<:Group, ω} = order(G)

FusionStyle(::Type{VecGωIrr{G, ω}}) where {G<:Group, ω}  = SimpleFusion()
BraidingStyle(::Type{VecGωIrr{G, ω}}) where {G<:Group, ω}  = NoBraiding()
Nsymbol(a::VecGωIrr{G, ω}, b::VecGωIrr{G, ω}, c::VecGωIrr{G, ω}) where {G<:Group, ω} = (c.g == a.g * b.g)
function Fsymbol(a::VecGωIrr{G, ω}, b::VecGωIrr{G, ω}, c::VecGωIrr{G, ω}, d::VecGωIrr{G, ω}, e::VecGωIrr{G, ω}, f::VecGωIrr{G, ω}) where {G<:Group, ω}
    return Nsymbol(a, b, e) * Nsymbol(e, c, d) * Nsymbol(b, c, f) * Nsymbol(a, f, d) * ω(a.g, b.g, c.g)
end
function Base.one(::Type{VecGωIrr{G, ω}}) where {G<:Group, ω}
    return VecGωIrr{G, ω}(identity_element(G))
end
Base.conj(c::VecGωIrr{G, ω}) where {G<:Group, ω} = VecGωIrr{G, ω}(inverse(c.g))
⊗(c1::VecGωIrr{G, ω}, c2::VecGωIrr{G, ω}) where {G<:Group, ω} = (VecGωIrr{G, ω}(c1.g*c2.g),)


Base.IteratorSize(::Type{SectorValues{VecGωIrr{G, ω}}}) where {G<:Group, ω} = rank(VecGωIrr{G, ω})==Inf ? IsInfinite() : HasLength()
Base.length(::SectorValues{VecGωIrr{G, ω}}) where {G<:Group, ω} = rank(VecGωIrr{G, ω})
Base.getindex(::SectorValues{VecGωIrr{G, ω}}, i::Int) where {G<:Group, ω} = VecGωIrr{G, ω}(G[i])

function Base.iterate(::SectorValues{VecGωIrr{G, ω}}, i::Int=0)  where {G<:Group, ω}
    if rank(VecGωIrr{G, ω})==Inf
        return i <= 0 ? (VecGωIrr{G, ω}(G[i]), (-i + 1)) : (VecGωIrr{G, ω}(G[i]), -i)
    else
        return i == rank(VecGωIrr{G, ω})-1 ? nothing : (VecGωIrr{G, ω}(G[i+1]), i + 1)
    end
end
findindex(::SectorValues{VecGωIrr{G, ω}}, g::VecGωIrr{G, ω})  where {G<:Group, ω} = findindex(g.g)

Base.isless(c1::VecGωIrr{G, ω}, c2::VecGωIrr{G, ω}) where {G<:Group, ω} = isless(findindex(c1.g), findindex(c2.g))


