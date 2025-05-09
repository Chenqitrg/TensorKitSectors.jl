# PlanarTrivial
"""
    struct PlanarTrivial <: Sector
    PlanarTrivial()

Represents a trivial anyon sector, i.e. a trivial sector without braiding. This is mostly
useful for testing.
"""
struct PlanarTrivial <: Sector end

Base.IteratorSize(::Type{SectorValues{PlanarTrivial}}) = HasLength()
Base.length(::SectorValues{PlanarTrivial}) = 1
Base.iterate(::SectorValues{PlanarTrivial}, i=0) = i == 0 ? (PlanarTrivial(), 1) : nothing
function Base.getindex(::SectorValues{PlanarTrivial}, i::Int)
    return i == 1 ? PlanarTrivial() :
           throw(BoundsError(values(PlanarTrivial), i))
end
findindex(::SectorValues{PlanarTrivial}, c::PlanarTrivial) = 1
Base.isless(::PlanarTrivial, ::PlanarTrivial) = false

Base.one(::Type{PlanarTrivial}) = PlanarTrivial()
Base.conj(::PlanarTrivial) = PlanarTrivial()

FusionStyle(::Type{PlanarTrivial}) = UniqueFusion()
BraidingStyle(::Type{PlanarTrivial}) = NoBraiding()
Base.isreal(::Type{PlanarTrivial}) = true

Nsymbol(::Vararg{PlanarTrivial,3}) = 1
Fsymbol(::Vararg{PlanarTrivial,6}) = 1

⊗(::PlanarTrivial, ::PlanarTrivial) = (PlanarTrivial(),)

# FibonacciAnyons
"""
    struct FibonacciAnyon <: Sector
    FibonacciAnyon(s::Symbol)

Represents the anyons of the Fibonacci modular fusion category. It can take two values,
corresponding to the trivial sector `FibonacciAnyon(:I)` and the non-trivial sector
`FibonacciAnyon(:τ)` with fusion rules ``τ ⊗ τ = 1 ⊕ τ``.

## Fields
- `isone::Bool`: indicates whether the sector corresponds the to trivial anyon `:I`
  (`true`), or the non-trivial anyon `:τ` (`false`).
"""
struct FibonacciAnyon <: Sector
    isone::Bool
    function FibonacciAnyon(s::Symbol)
        s in (:I, :τ, :tau) || throw(ArgumentError("Unknown FibonacciAnyon $s."))
        return new(s === :I)
    end
end

Base.IteratorSize(::Type{SectorValues{FibonacciAnyon}}) = HasLength()
Base.length(::SectorValues{FibonacciAnyon}) = 2
function Base.iterate(::SectorValues{FibonacciAnyon}, i=0)
    return i == 0 ? (FibonacciAnyon(:I), 1) : (i == 1 ? (FibonacciAnyon(:τ), 2) : nothing)
end
function Base.getindex(S::SectorValues{FibonacciAnyon}, i::Int)
    if i == 1
        return FibonacciAnyon(:I)
    elseif i == 2
        return FibonacciAnyon(:τ)
    else
        throw(BoundsError(S, i))
    end
end
findindex(::SectorValues{FibonacciAnyon}, s::FibonacciAnyon) = 2 - s.isone

Base.convert(::Type{FibonacciAnyon}, s::Symbol) = FibonacciAnyon(s)
Base.one(::Type{FibonacciAnyon}) = FibonacciAnyon(:I)
Base.conj(s::FibonacciAnyon) = s

const _goldenratio = Float64(MathConstants.golden)
dim(a::FibonacciAnyon) = isone(a) ? one(_goldenratio) : _goldenratio

FusionStyle(::Type{FibonacciAnyon}) = SimpleFusion()
BraidingStyle(::Type{FibonacciAnyon}) = Anyonic()
Base.isreal(::Type{FibonacciAnyon}) = false

⊗(a::FibonacciAnyon, b::FibonacciAnyon) = FibonacciIterator(a, b)

struct FibonacciIterator
    a::FibonacciAnyon
    b::FibonacciAnyon
end
Base.IteratorSize(::Type{FibonacciIterator}) = Base.HasLength()
Base.IteratorEltype(::Type{FibonacciIterator}) = Base.HasEltype()
Base.length(iter::FibonacciIterator) = (isone(iter.a) || isone(iter.b)) ? 1 : 2
Base.eltype(::Type{FibonacciIterator}) = FibonacciAnyon
function Base.iterate(iter::FibonacciIterator, state=1)
    I = FibonacciAnyon(:I)
    τ = FibonacciAnyon(:τ)
    if state == 1 # first iteration
        iter.a == I && return (iter.b, 2)
        iter.b == I && return (iter.a, 2)
        return (I, 2)
    elseif state == 2
        (iter.a == iter.b == τ) && return (τ, 3)
        return nothing
    else
        return nothing
    end
end

function Nsymbol(a::FibonacciAnyon, b::FibonacciAnyon, c::FibonacciAnyon)
    return isone(a) + isone(b) + isone(c) != 2
end # zero if one tau and two ones

function Fsymbol(a::FibonacciAnyon, b::FibonacciAnyon, c::FibonacciAnyon,
                 d::FibonacciAnyon, e::FibonacciAnyon, f::FibonacciAnyon)
    Nsymbol(a, b, e) || return zero(_goldenratio)
    Nsymbol(e, c, d) || return zero(_goldenratio)
    Nsymbol(b, c, f) || return zero(_goldenratio)
    Nsymbol(a, f, d) || return zero(_goldenratio)

    I = FibonacciAnyon(:I)
    τ = FibonacciAnyon(:τ)
    if a == b == c == d == τ
        if e == f == I
            return +1 / _goldenratio
        elseif e == f == τ
            return -1 / _goldenratio
        else
            return +1 / sqrt(_goldenratio)
        end
    else
        return one(_goldenratio)
    end
end

function Rsymbol(a::FibonacciAnyon, b::FibonacciAnyon, c::FibonacciAnyon)
    Nsymbol(a, b, c) || return 0 * cispi(0 / 1)
    if isone(a) || isone(b)
        return cispi(0 / 1)
    else
        return isone(c) ? cispi(4 / 5) : cispi(-3 / 5)
    end
end

function Base.show(io::IO, a::FibonacciAnyon)
    s = isone(a) ? ":I" : ":τ"
    return get(io, :typeinfo, nothing) === FibonacciAnyon ?
           print(io, s) : print(io, "FibonacciAnyon(", s, ")")
end

Base.hash(a::FibonacciAnyon, h::UInt) = hash(a.isone, h)
Base.isless(a::FibonacciAnyon, b::FibonacciAnyon) = isless(!a.isone, !b.isone)

# IsingAnyons
"""
    struct IsingAnyon <: Sector
    IsingAnyon(s::Symbol)

Represents the anyons of the Ising modular fusion category. It can take three values,
corresponding to the trivial sector `IsingAnyon(:I)` and the non-trivial sectors
`IsingAnyon(:σ)` and `IsingAnyon(:ψ)`, with fusion rules ``ψ ⊗ ψ = 1``, ``σ ⊗ ψ = σ``, and
``σ ⊗ σ = 1 ⊕ ψ``.

## Fields
- `s::Symbol`: the label of the represented anyon, which can be `:I`, `:σ`, or `:ψ`.
"""
struct IsingAnyon <: Sector
    s::Symbol
    function IsingAnyon(s::Symbol)
        s == :sigma && (s = :σ)
        s == :psi && (s = :ψ)
        if !(s in (:I, :σ, :ψ))
            throw(ArgumentError("Unknown IsingAnyon $s."))
        end
        return new(s)
    end
end

const all_isinganyons = (IsingAnyon(:I), IsingAnyon(:σ), IsingAnyon(:ψ))

Base.IteratorSize(::Type{SectorValues{IsingAnyon}}) = HasLength()
Base.length(::SectorValues{IsingAnyon}) = length(all_isinganyons)
Base.iterate(::SectorValues{IsingAnyon}, i=1) = iterate(all_isinganyons, i)
Base.getindex(::SectorValues{IsingAnyon}, i::Int) = getindex(all_isinganyons, i)

function findindex(::SectorValues{IsingAnyon}, a::IsingAnyon)
    a == all_isinganyons[1] && return 1
    a == all_isinganyons[2] && return 2
    return 3
end

Base.convert(::Type{IsingAnyon}, s::Symbol) = IsingAnyon(s)
Base.one(::Type{IsingAnyon}) = IsingAnyon(:I)
Base.conj(s::IsingAnyon) = s

dim(a::IsingAnyon) = a.s == :σ ? sqrt(2) : 1.0

FusionStyle(::Type{IsingAnyon}) = SimpleFusion()
BraidingStyle(::Type{IsingAnyon}) = Anyonic()
Base.isreal(::Type{IsingAnyon}) = false

⊗(a::IsingAnyon, b::IsingAnyon) = IsingIterator(a, b)

struct IsingIterator
    a::IsingAnyon
    b::IsingAnyon
end

Base.IteratorSize(::Type{IsingIterator}) = Base.HasLength()
Base.IteratorEltype(::Type{IsingIterator}) = Base.HasEltype()
Base.eltype(::Type{IsingIterator}) = IsingAnyon

function Base.length(iter::IsingIterator)
    σ = IsingAnyon(:σ)
    return (iter.a == σ && iter.b == σ) ? 2 : 1
end

function Base.iterate(iter::IsingIterator, state=1)
    I, σ, ψ = all_isinganyons
    if state == 1 # first iteration
        iter.a == I && return (iter.b, 2)
        iter.b == I && return (iter.a, 2)
        iter.a == σ && iter.b == ψ && return (σ, 2)
        iter.a == ψ && iter.b == σ && return (σ, 2)
        return (I, 2)
    elseif state == 2
        (iter.a == iter.b == σ) && return (ψ, 3)
        return nothing
    else
        return nothing
    end
end

function Nsymbol(a::IsingAnyon, b::IsingAnyon, c::IsingAnyon)
    I, σ, ψ = all_isinganyons
    return ((a == I && b == c)
            || (b == I && a == c)
            || (c == I && a == b)
            || (a == σ && b == σ && c == ψ)
            || (a == σ && b == ψ && c == σ)
            || (a == ψ && b == σ && c == σ))
end

function Fsymbol(a::IsingAnyon, b::IsingAnyon, c::IsingAnyon,
                 d::IsingAnyon, e::IsingAnyon, f::IsingAnyon)
    Nsymbol(a, b, e) || return 0.0
    Nsymbol(e, c, d) || return 0.0
    Nsymbol(b, c, f) || return 0.0
    Nsymbol(a, f, d) || return 0.0
    I, σ, ψ = all_isinganyons
    if a == b == c == d == σ
        if e == f == ψ
            return -1.0 / sqrt(2.0)
        else
            return 1.0 / sqrt(2.0)
        end
    end
    if e == f == σ
        if a == c == σ && b == d == ψ
            return -1.0
        elseif a == c == ψ && b == d == σ
            return -1.0
        end
    end
    return 1.0
end

function Rsymbol(a::IsingAnyon, b::IsingAnyon, c::IsingAnyon)
    Nsymbol(a, b, c) || return complex(0.0)
    I, σ, ψ = all_isinganyons
    if c == I
        if b == a == σ
            return cispi(-1 / 8)
        elseif b == a == ψ
            return complex(-1.0)
        end
    elseif c == σ && (a == σ && b == ψ || a == ψ && b == σ)
        return -1.0im
    elseif c == ψ && a == b == σ
        return cispi(3 / 8)
    end
    return complex(1.0)
end

function Base.show(io::IO, a::IsingAnyon)
    if get(io, :typeinfo, nothing) === IsingAnyon
        return print(io, ":$(a.s)")
    else
        return print(io, "IsingAnyon(:$(a.s))")
    end
end

Base.hash(s::IsingAnyon, h::UInt) = hash(s.s, h)

function Base.isless(a::IsingAnyon, b::IsingAnyon)
    vals = SectorValues{IsingAnyon}()
    return isless(findindex(vals, a), findindex(vals, b))
end



"""
struct Semion <: Sector
    isone::Bool
    function Semion(s::Integer)
        s in (0, 1) || throw(ArgumentError("Unknown Semion s."))
        return new(s === 0)
    end
end
The semion category is almost Rep(ℤ₂), except that its associator is non-trivial and has non-trivial braiding. 
- The associator is non-trivial only when three fusion anyons are all semions: ω(1,1,1) = -1.
- And the braiding is non-trivial only when two fusion anyons are semions: R(1,1) = i.
Moreover, it is, perhaps the simplest, non-trivial unitary modular tensor category. It has central charge c = 1.
It has two simple objects, the trivial one, denoted by 0 and the semion, denoted by 1, which is a non-trivial anyon.

Moreover, the Frobenius-Schur indicator is -1 for the semion and +1 for the trivial anyon.
"""
struct Semion <: Sector
    isone::Bool
    function Semion(s::Integer)
        s in (0, 1) || throw(ArgumentError("Unknown Semion $s."))
        return new(s === 0)
    end
end

Base.IteratorSize(::Type{SectorValues{Semion}}) = Base.HasLength()
Base.length(::SectorValues{Semion}) = 2
function Base.iterate(::SectorValues{Semion}, i=0)
    return i == 0 ? (Semion(0), 1) : (i == 1 ? (Semion(1), 2) : nothing)
end
function Base.getindex(::SectorValues{Semion}, i::Int)
    return 1 <= i <= 2 ? Semion(i - 1) : throw(BoundsError(values(Semion), i))
end
findindex(::SectorValues{Semion}, s::Semion) = 2 - s.isone

Base.convert(::Type{Semion}, s::Int) = Semion(s)
Base.one(::Type{Semion}) = Semion(0)
Base.conj(s::Semion) = s
frobeniusschur(s::Semion) = s==Semion(0) ? 1 : -1

dim(a::Semion) = 1

FusionStyle(::Type{Semion}) = SimpleFusion()
BraidingStyle(::Type{Semion}) = Anyonic()
Base.isreal(::Type{Semion}) = false

⊗(a::Semion, b::Semion) = (Semion(mod(to_Z2(a) + to_Z2(b),2)),)


Base.IteratorSize(::SectorValues{Semion}) = Base.HasLength()
Base.IteratorEltype(::SectorValues{Semion}) = Base.HasEltype()
Base.length(::SectorValues{Semion}) = 2
function Base.iterate(::SectorValues{Semion}, i=0)
    return i == 2 ? nothing : (Semion(i), i + 1)
end


function to_Z2(s::Semion)
    return s.isone ? 0 : 1
end
function Nsymbol(a::Semion, b::Semion, c::Semion)
    return mod(to_Z2(a) + to_Z2(b) + to_Z2(c),2) == 0 ? 1 : 0
end

function Fsymbol(a::Semion, b::Semion, c::Semion,
                 d::Semion, e::Semion, f::Semion)
    if a == b == c == d == Semion(1)
        return Int(- Nsymbol(a, b, e) * Nsymbol(e, c, d) * Nsymbol(b, c, f) * Nsymbol(a, f, d))
    else
        return Int(Nsymbol(a, b, e) * Nsymbol(e, c, d) * Nsymbol(b, c, f) * Nsymbol(a, f, d))
    end
end

function Rsymbol(a::Semion, b::Semion, c::Semion)
    if isone(a) && isone(b)
        return Complex(im* Nsymbol(a, b, c))
    else
        return Complex(Nsymbol(a, b, c))
    end
end

function Base.show(io::IO, a::Semion)
    s = isone(a) ? "0" : "1"
    return get(io, :typeinfo, nothing) === Semion ?
           print(io, s) : print(io, "Semion(", s, ")")
end

Base.hash(a::Semion, h::UInt) = hash(a.isone, h)
Base.isless(a::Semion, b::Semion) = isless(!a.isone, !b.isone)


