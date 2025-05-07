struct Irr{𝒞<:FusionCategory} <: Sector
    value::Any
    function Irr{𝒞}(value) where {𝒞<:VecGω}
        G = 𝒞.parameters[1]
        if value isa GroupElement{G}
            new(value)
        else
            throw(ArgumentError("Irr value must be a GroupElement of type $G"))
        end
    end
end

rank(::Type{𝒞}) where {𝒞<:VecGω} = order(𝒞.parameters[1])

FusionStyle(::Type{Irr{𝒞}}) where {𝒞<:VecGω}  = UniqueFusion()
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
findindex(::SectorValues{Irr{𝒞}}, g::Irr{𝒞})  where {𝒞<:VecGω} = findindex(𝒞.parameters[1], g.value)

f(a, b, c) = 1
# FusionStyle(Irr{VecGω{ℤ{3}, f}})

rank(VecGω{ℤ{3}×D{5}, f})

ℤ{Inf}[-5]
(ℤ{Inf}×ℤ{5})[-15]
findindex(D{3}, GroupElement{D{3}}(0,0))
# @show GroupElement{ℤ{3}}(0)
# @show CohomologyGroup{3, ℤ{3}, ℤ{Inf}}
# GroupElement{CohomologyGroup{3, ℤ{3}, ℤ{Inf}}}(f)
# Irr{VecGω{ℤ{3}, ω}}(GroupElement{ℤ{3}}(1))


