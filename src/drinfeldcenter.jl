using TensorKit
using TensorKitSectors

abstract type ℨ{𝒞<:Sector} <: Sector end
struct QDZ{N} <: ℨ{ZNIrrep{N}}
    charge::Int
    flux::Int
    function QDZ{N}(e::Int, m::Int) where {N}
        new{N}(mod(e,N), mod(m, N))
    end
end 
FusionStyle(::Type{QDZ{N}}) where {N}  = SimpleFusion()
BraidingStyle(::Type{QDZ{N}}) where {N}  = Anyonic()
Nsymbol(a::QDZ{N}, b::QDZ{N}, c::QDZ{N}) where {N} = (mod(a.charge+b.charge, N)==c.charge) && (mod(a.flux+b.flux, N)==c.flux)
Fsymbol(a::QDZ{N}, b::QDZ{N}, c::QDZ{N}, d::QDZ{N}, e::QDZ{N}, f::QDZ{N}) where {N} = Nsymbol(a, b, e) * Nsymbol(e, c, d) * Nsymbol(b, c, f) * Nsymbol(a, f, d)
Base.one(::Type{QDZ{N}}) where {N} = QDZ{N}(0,0)
Base.conj(c::QDZ{N}) where {N} = QDZ{N}(-c.charge, -c.flux)
⊗(c1::QDZ{N}, c2::QDZ{N}) where {N} = (QDZ{N}(c1.charge+c2.charge, c1.flux+c2.flux),)
Base.IteratorSize(::Type{SectorValues{QDZ{N}}}) where {N} = HasLength()
Base.length(::SectorValues{QDZ{N}}) where {N} = N^2
Base.getindex(::SectorValues{QDZ{N}}, e::Int, m::Int) where {N} = QDZ{N}(e-1,m-1)
function Base.iterate(::SectorValues{QDZ{N}}, i::Int=1)  where {N}
    if rank(VecGωIrr{G, ω})==Inf
        return i <= 1 ? (SectorValues{VecGωIrr{G, ω}}()[i], (-i + 1)) : (VecGωIrr{G, ω}[i], -i)
    else
        return i == rank(VecGωIrr{G, ω}) ? nothing : (VecGωIrr{G, ω}[i], i + 1)
    end
end

# SectorValues{QDZ{3}}()[1,2]