# Groups
#------------------------------------------------------------------------------#
abstract type Group end
abstract type AbelianGroup <: Group end

abstract type ℤ{N} <: AbelianGroup end
struct D{N} <: Group end
abstract type U₁ <: AbelianGroup end
abstract type SU{N} <: Group end
abstract type CU₁ <: Group end

const ℤ₂ = ℤ{2}
const ℤ₃ = ℤ{3}
const ℤ₄ = ℤ{4}
const SU₂ = SU{2}
const D₃ = D{3}
const D₄ = D{4}

type_repr(::Type{ℤ₂}) = "ℤ₂"
type_repr(::Type{ℤ₃}) = "ℤ₃"
type_repr(::Type{ℤ₄}) = "ℤ₄"
type_repr(::Type{ℤ{Inf}}) = "ℤ"
type_repr(::Type{SU₂}) = "SU₂"
type_repr(T::Type) = repr(T)

const GroupTuple = Tuple{Vararg{Group}}

abstract type ProductGroup{T<:GroupTuple} <: Group end

"""
    ×(G::Vararg{Type{<:Group}}) -> ProductGroup{Tuple{G...}}
    times(G::Vararg{Type{<:Group}}) -> ProductGroup{Tuple{G...}}

Construct the direct product of a (list of) groups.
"""
function ×(::Vararg{Type{<:Group}}) end
const times = ×

×(a::Type{<:Group}, b::Type{<:Group}, c::Type{<:Group}...) = ×(×(a, b), c...)
×(G::Type{<:Group}) = ProductGroup{Tuple{G}}
×(G1::Type{ProductGroup{Tuple{}}},
G2::Type{ProductGroup{T}}) where {T<:GroupTuple} = G2
function ×(G1::Type{ProductGroup{T1}},
           G2::Type{ProductGroup{T2}}) where {T1<:GroupTuple,T2<:GroupTuple}
    return tuple_type_head(T1) × (ProductGroup{tuple_type_tail(T1)} × G2)
end
×(G1::Type{ProductGroup{Tuple{}}}, G2::Type{<:Group}) = ProductGroup{Tuple{G2}}
function ×(G1::Type{ProductGroup{T}}, G2::Type{<:Group}) where {T<:GroupTuple}
    return Base.tuple_type_head(T) × (ProductGroup{Base.tuple_type_tail(T)} × G2)
end
function ×(G1::Type{<:Group}, G2::Type{ProductGroup{T}}) where {T<:GroupTuple}
    return ProductGroup{Base.tuple_type_cons(G1, T)}
end
×(G1::Type{<:Group}, G2::Type{<:Group}) = ProductGroup{Tuple{G1,G2}}


order(::Type{ℤ{N}}) where {N} = N
order(::Type{D{N}}) where {N} = 2*N
function order(::Type{ProductGroup{Gs}}) where {Gs<:GroupTuple}
    orders = map(order, Gs.parameters)
    return prod(orders)
end

abstract type FusionCategory end
abstract type VecGω{G,ω} <: FusionCategory end
abstract type CohomologyGroup{N,G,A} <: AbelianGroup end
abstract type TY{A, χ, ϵ} <: FusionCategory end # A: an Abelian group, χ: a symmetric non-degenerate bi-character, ϵ: the Frobenius-Schur indicator for the non-invertible object :σ

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
    function Irr{𝒞}(obj) where {𝒞<:TY}
        A = 𝒞.parameters[1]
        if isa(obj, GroupElement{A}) || obj == :σ
            new{𝒞}(obj)
        else
            throw(ArgumentError("Illegal object $obj"))
        end
    end
end


