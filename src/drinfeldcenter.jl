abstract type ℨ{𝒞<:Sector} <: ModularSector end

"""
    struct QDℤ{N} <: ℨ{ZNIrrep{N}}
        charge::Int
        flux::Int
    end
    
The quantum double of Rep(ℤₙ) or Vec(ℤₙ).
"""
struct QDℤ{N} <: ℨ{ZNIrrep{N}}
    charge::Int
    flux::Int
    function QDℤ{N}(e::Int, m::Int) where {N}
        new{N}(mod(e,N), mod(m, N))
    end
end

take_center(::Type{ZNIrrep{N}}) where {N} = QDℤ{N}
sector_rank(::Type{QDℤ{N}}) where {N} = N^2

FusionStyle(::Type{QDℤ{N}}) where {N}  = SimpleFusion()
BraidingStyle(::Type{QDℤ{N}}) where {N}  = Anyonic()
is_modular(::Type{QDℤ{N}}) where {N} = true

Nsymbol(a::QDℤ{N}, b::QDℤ{N}, c::QDℤ{N}) where {N} = (mod(a.charge+b.charge, N)==c.charge) && (mod(a.flux+b.flux, N)==c.flux)
Fsymbol(a::QDℤ{N}, b::QDℤ{N}, c::QDℤ{N}, d::QDℤ{N}, e::QDℤ{N}, f::QDℤ{N}) where {N} = Nsymbol(a, b, e) * Nsymbol(e, c, d) * Nsymbol(b, c, f) * Nsymbol(a, f, d)
Base.one(::Type{QDℤ{N}}) where {N} = QDℤ{N}(0,0)
Base.conj(c::QDℤ{N}) where {N} = QDℤ{N}(-c.charge, -c.flux)
⊗(c1::QDℤ{N}, c2::QDℤ{N}) where {N} = (QDℤ{N}(c1.charge+c2.charge, c1.flux+c2.flux),)
Base.IteratorSize(::Type{SectorValues{QDℤ{N}}}) where {N} = HasLength()
Base.length(::SectorValues{QDℤ{N}}) where {N} = N^2
function Base.getindex(::SectorValues{QDℤ{N}}, i::Int) where {N}
    e = (i-1)÷N
    m = (i-1)%N
    return QDℤ{N}(e, m)
end
Base.iterate(::SectorValues{QDℤ{N}}, i::Int=1)  where {N} = i == N^2+1 ? nothing : (QDℤ{N}((i-1)÷N,(i-1)%N), i + 1)
findindex(::SectorValues{QDℤ{N}}, a::QDℤ{N})  where {N} = a.charge*N + a.flux + 1
function Base.isless(c1::QDℤ{N}, c2::QDℤ{N}) where {N}
    if c1.charge < c2.charge
        return true
    elseif (c1.charge == c2.charge) && c1.flux < c2.flux
        return true
    else
        return false
    end
end
function Rsymbol(a::QDℤ{N}, b::QDℤ{N}, c::QDℤ{N}) where {N}
    R = Nsymbol(a, b, c) * exp(2 * pi * im/N * a.flux * b.charge)
    return R
end


"""
    struct QDAb{A} <: ℨ{VecGIrr{A}}
        a::A
        f::χ{A}
    end

The quantum double of Vec_A, where A is an Abelian group.

Its simple objects are elements in A and elements in A*.

"""
struct QDAb{A} <: ℨ{VecGIrr{A}}
    a::A
    f::χ{A}
    function QDAb{A}(a::A, f::χ{A}) where {A<:Group}
        if !is_abelian(A)
            throw(ArgumentError("The group $A is not Abelian"))
        end
        new{A}(a,f)
    end
end

take_center(::Type{VecGIrr{A}}) where {A<:Group} = QDAb{A}

FusionStyle(::Type{QDAb{A}}) where {A<:Group}  = SimpleFusion()
BraidingStyle(::Type{QDAb{A}}) where {A<:Group}  = Anyonic()
is_modular(::Type{QDAb{A}}) where {A<:Group}  = true

function Nsymbol(a::QDAb{A}, b::QDAb{A}, c::QDAb{A}) where {A<:Group}
    if a.a*b.a != c.a
        return false
    elseif a.f * b.f != c.f
        return false
    end
    return true
end
Fsymbol(a::QDAb{A}, b::QDAb{A}, c::QDAb{A}, d::QDAb{A}, e::QDAb{A}, f::QDAb{A}) where {A<:Group} = Nsymbol(a, b, e) * Nsymbol(e, c, d) * Nsymbol(b, c, f) * Nsymbol(a, f, d)
Base.one(::Type{QDAb{A}}) where {A<:Group} = QDAb{A}(identity_element(A),identity_element(χ{A}))
Base.conj(c::QDAb{A}) where {A<:Group} = QDAb{A}(inverse(c.a), inverse(c.f))
⊗(c1::QDAb{A}, c2::QDAb{A}) where {A<:Group} = (QDAb{A}(c1.a*c2.a, c1.f*c2.f),)
function Rsymbol(a::QDAb{A}, b::QDAb{A}, c::QDAb{A}) where {A<:Group}
    R = Nsymbol(a, b, c) * eval(a.f, b.a)
    return R
end

# TO DO: choose a basis for character such that the iterator can be defined

forget_flux(a::QDAb{A}) where {A<:Group} = VecGIrr{A}(a.a)


"""
The Drinfeld center of Tambara-Yamagami category is given in the paper arXiv:0905.3117
"""
struct ℨTY{A,χ,ϵ}<:ℨ{TYIrr{A,χ,ϵ}}

end

take_center(::Type{𝒞}) where {𝒞<:ModularSector} = 𝒞 ⊠ TimeReversed{𝒞}