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
function Base.getindex(::SectorValues{QDZ{N}}, i::Int) where {N}
    e = (i-1)÷N
    m = (i-1)%N
    return QDZ{N}(e, m)
end
Base.iterate(::SectorValues{QDZ{N}}, i::Int=1)  where {N} = i == N^2 ? nothing : (QDZ{N}((i-1)÷N,(i-1)%N), i + 1)
findindex(::SectorValues{QDZ{N}}, a::QDZ{N})  where {N} = a.charge*N + a.flux + 1
function Base.isless(c1::QDZ{N}, c2::QDZ{N}) where {N}
    if c1.charge < c2.charge
        return true
    elseif (c1.charge == c2.charge) && c1.flux < c2.flux
        return true
    else
        return false
    end
end

function Rsymbol(a::QDZ{N}, b::QDZ{N}, c::QDZ{N}) where {N}
    R = Nsymbol(a, b, c) * exp(2 * pi * im/N * a.flux * b.charge)
    return R
end

forget_flux(a::QDZ{N}) where {N} = ZNIrrep{N}(a.charge)

foget_charge(a::QDZ{N}) where {N} = ZNIrrep{N}(a.flux)

function HalfBraiding_charge(a::QDZ{N}, V::GradedSpace{ZNIrrep{N}, NTuple{N, Int64}}) where {N}
    fgt_a = forget_flux(a)
    W = ZNSpace{N}(fgt_a=>1)
    Ω = zeros(V⊗W←W⊗V)
    for i = 0:N-1
        
    end
end
# iden = QDZ{2}(0,0)
# e = QDZ{2}(1,0)
# m = QDZ{2}(0,1)
# f = QDZ{2}(1,1)
# Vect[QDZ{2}](e=>1,m=>1)