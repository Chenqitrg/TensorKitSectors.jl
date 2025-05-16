# Groups
#------------------------------------------------------------------------------#
abstract type Group end
abstract type AbelianGroup <: Group end

struct ℤ{N} <: AbelianGroup 
    a::Int
    function ℤ{N}(n::Int) where {N}
        if N == Inf
            throw("The group order $N must be finite")
        end

        new{N}(mod(n, N))
    end
end
struct D{N} <: Group 
    s::Int
    r::Int
    function D{N}(s::Int, r::Int) where {N}
        if N == Inf
            throw("The group order 2×$N must be finite")
        end
        new{N}(mod(s, 2), mod(r, N))
    end
end

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
type_repr(::Type{SU₂}) = "SU₂"
type_repr(T::Type) = repr(T)

const GroupTuple = Tuple{Vararg{Group}}

struct ProductGroup{T} <: Group 
    components
    function ProductGroup{T}(elements...) where {T<:GroupTuple}
        for (g, G) in zip(elements, T.parameters)
            if !(g isa G)
                throw(ArgumentError("The element $g is not in group $G"))
            end
        end
        new{T}(elements)
    end
end

"""
    ×(G::Vararg{Type{<:Group}}) -> ProductGroup{Tuple{G...}}
    times(G::Vararg{Type{<:Group}}) -> ProductGroup{Tuple{G...}}

Construct the direct product of a (list of) groups.
"""
function ×(::Vararg{Type{<:Group}}) end # Becareful that the LinearAlgebra also exports "×"
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

# Default trait: assume any group is non-abelian.
_is_abelian(::Type{G}) where {G<:Group} = false

# Mark all subtypes of AbelianGroup as abelian.
_is_abelian(::Type{G}) where {G<:AbelianGroup} = true

# For a ProductGroup, check if all component groups are abelian.
function _is_abelian(::Type{ProductGroup{T}}) where {T<:GroupTuple}
    # T.parameters returns a tuple of the component types.
    # We use `all` to verify each component is abelian.
    all(g -> _is_abelian(g), T.parameters)
end

# Public function to test if a group type is abelian.
is_abelian(::Type{G}) where {G<:Group} = _is_abelian(G)


# ===========================================
# Group Operations
# =========================================== 

function elements(::Type{ℤ{N}}) where {N}
    return ntuple(i -> ℤ{N}(i - 1), N)
end

# elements(ℤ{3})
function elements(::Type{D{N}}) where {N}
    rotations = ntuple(i -> D{N}(0, i - 1), N)  # (e, r, r², ...)
    reflections = ntuple(i -> D{N}(1, i - 1), N)  # (s, sr, sr², ...)
    return (rotations..., reflections...)  # Joining the two tuples
end

# elements(D{3})
function elements(::Type{ProductGroup{Gs}}) where {Gs<:GroupTuple}
    group_elements = map(elements, Gs.parameters)
    cartesian_product = collect(Iterators.product(group_elements...))
    return tuple(map(x->ProductGroup{Gs}(x...),cartesian_product)...)
end

# display(elements(ℤ{3}×D{3}))

identity_element(::Type{ℤ{N}}) where {N} = ℤ{N}(0)
identity_element(::Type{D{N}}) where {N} = D{N}(0,0)
function identity_element(::Type{ProductGroup{Gs}}) where {Gs<:GroupTuple}
    groups = Gs.parameters
    return ProductGroup{Gs}((identity_element(G) for G in groups)...)
end

function inverse(x::ℤ{N}) where {N}
    return ℤ{N}(-x.a)
end
function inverse(x::D{N}) where {N}
    return D{N}(-x.s, (-1)^(x.s + 1) * x.r)
end
function inverse(x::ProductGroup{Gs}) where {Gs<:GroupTuple}
    groups = Gs.parameters[1].parameters
    inverse_elements = map(p -> inverse(p), x.components)
    return ProductGroup{Gs}(inverse_elements)
end

function Base.:*(x::ℤ{N}, y::ℤ{N}) where {N}
    return ℤ{N}(mod(x.a + y.a, N))
end
function Base.:*(x::D{N}, y::D{N}) where {N}
    s1, r1 = x.s, x.r
    s2, r2 = y.s, y.r
    return D{N}(mod(s1 + s2, 2), mod((-1)^s2 * r1 + r2, N))
end
function Base.:*(x::ProductGroup{Gs}, y::ProductGroup{Gs}) where {Gs<:GroupTuple}
    newelement = ()
    for (x1, y1) in zip(x.components, y.components)
        newelement = (newelement..., x1 * y1)
    end
    return ProductGroup{Gs}(newelement...)
end

Base.getindex(::Type{ℤ{N}}, i::Int) where {N} = ℤ{N}(i-1) # Count from 1
Base.getindex(::Type{D{N}}, i::Int) where {N} = D{N}((i-1)÷N, (i-1)%N) # Count from 1
function Base.getindex(::Type{ProductGroup{Gs}}, i::Int) where {Gs<:GroupTuple}
    i -= 1
    groups = Gs.parameters
    elementtuple = ()
    for n in length(groups):-1:1
        group_tail = groups[n]
        order_tail = order(group_tail)
        ind_tail =  i%order_tail
        i = i÷order_tail
        elementtuple = (group_tail[ind_tail+1], elementtuple...)
    end
    return ProductGroup{Gs}(elementtuple...)
end

findindex(g::ℤ{N}) where {N} = g.a + 1 # Count from 1
findindex(g::D{N}) where {N} = g.s * N + g.r+1 # Count from 1
function findindex(g::ProductGroup{Gs}) where {Gs<:GroupTuple}
    index = 0
    groups = Gs.parameters
    weight = 1
    for i in length(groups):-1:1
        g_this = g.components[i]
        G_this = groups[i]
        index += (findindex(g_this)-1) * weight
        weight *= order(G_this)
    end
    return index+1
end

"""
    struct χ{A} <: AbelianGroup
        f
    end

The dual group of an Abelian group A.
"""
struct χ{A} <: AbelianGroup
    f
    function Fun{A}(f) where {A<:Group}
        if !is_abelian(A)
            throw(ArgumentError("$A is not Abelian"))
        end
        for a in elements(A), b in elements(A)
            if f(a*b) != f(a) * f(b)
                throw(ArgumentError("$f is not a character"))
            end
        end
        new{A}(f)
    end
end

identity_element(::Type{χ{A}}) where {A<:Group} = χ{A}(x::A -> 1.)
inverse(f::χ{A}) where {A<:Group} = χ{A}(x::A -> inv(f.f(x)))
function Base.:*(a::χ{A}, b::χ{A}) where {A<:Group}
    fab(g::A) = a.f(g) * b.f(g)
    return χ{A}(fab)
end

eval(f::χ{A}, a::A) where {A<:Group} = f.f(a)

import Base: ==
function ==(f::χ{A}, g::χ{A}) where {A<:Group} 
    for a in elements(A)
        if f.f(a) != g.f(a)
            return false
        end
    end
    return true
end

# Examples

# ℤ{3}(2)
# ℤ{3}(3)
# ℤ{3}(2) isa ℤ{3}
# ℤ{3}(1) * ℤ{3}(2)

# (ℤ₃×ℤ₃)(ℤ{3}(2), ℤ{3}(2))

# D₃(0,1) * D₃(1,0)

# inverse(D₃(0,1))

# identity_element(D₃)

# (ℤ₃×D₃)(ℤ₃(2), D₃(0,1))

# (D₃×ℤ₃)(D₃(0,1), ℤ₃(2)) * (D₃×ℤ₃)(D₃(1,1), ℤ₃(1))

# D₃[6]


# (D₃×ℤ₃)[7]

# findindex((D₃×ℤ₃)[7]) == 7

# findindex((D₃×ℤ₃)[4]) == 4

# order(ℤ₃)
# order(D₃×ℤ₃)

# G = ℤ₂×ℤ₂×ℤ₂
# x = ℤ₂(0)
# y = ℤ₂(1)
# G(y,y,y) * G(y,x,y) 


# elements(ℤ₃)
# elements(D₄)
# identity_element(ℤ₂×ℤ₂)
# elements(ℤ₂×ℤ₂)

# elements(D₃×ℤ₃)