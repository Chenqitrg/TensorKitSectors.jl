rank(::Type{𝒞}) where {𝒞<:VecGω} = order(𝒞.parameters[1])

FusionStyle(::Type{Irr{𝒞}}) where {𝒞<:VecGω}  = SimpleFusion()
BraidingStyle(::Type{Irr{𝒞}}) where {𝒞<:VecGω}  = NoBraiding()
Nsymbol(a::Irr{𝒞}, b::Irr{𝒞}, c::Irr{𝒞}) where {𝒞<:VecGω} = (c.value == a.value * b.value)
function Fsymbol(a::Irr{𝒞}, b::Irr{𝒞}, c::Irr{𝒞}, d::Irr{𝒞}, e::Irr{𝒞}, f::Irr{𝒞}) where {𝒞<:VecGω}
    ω = 𝒞.parameters[2]
    return Nsymbol(a, b, e) * Nsymbol(e, c, d) * Nsymbol(b, c, f) * Nsymbol(a, f, d) * ω(a.value, b.value, c.value)
end
function Base.one(::Type{Irr{𝒞}}) where {𝒞<:VecGω}
    G = 𝒞.parameters[1]
    return Irr{𝒞}(identity_element(G))
end
Base.conj(c::Irr{𝒞}) where {𝒞<:VecGω} = Irr{𝒞}(inverse(c.value))
⊗(c1::Irr{𝒞}, c2::Irr{𝒞}) where {𝒞<:VecGω} = (Irr{𝒞}(c1.value*c2.value),)


Base.IteratorSize(::Type{SectorValues{Irr{𝒞}}}) where {𝒞<:VecGω} = rank(𝒞)==Inf ? IsInfinite() : HasLength()
Base.length(::SectorValues{Irr{𝒞}}) where {𝒞<:VecGω} = rank(𝒞)
Base.getindex(::SectorValues{Irr{𝒞}}, i::Int) where {𝒞<:VecGω} = Irr{𝒞}(𝒞.parameters[1][i])

function Base.iterate(::SectorValues{Irr{𝒞}}, i::Int=0)  where {𝒞<:VecGω}
    if rank(𝒞)==Inf
        return i <= 0 ? (Irr{𝒞}[i], (-i + 1)) : (Irr{𝒞}[i], -i)
    else
        return i == rank(𝒞) ? nothing : (Irr{𝒞}[i], i + 1)
    end
end
findindex(::SectorValues{Irr{𝒞}}, g::Irr{𝒞})  where {𝒞<:VecGω} = findindex(g.value)


# Base.hash(c::ZNIrrep{N}, h::UInt) where {N} = hash(c.n, h)
Base.isless(c1::Irr{𝒞}, c2::Irr{𝒞}) where {𝒞<:VecGω} = isless(findindex(c1.value), findindex(c2.value))


