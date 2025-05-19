struct TimeReversed{I<:BraidedSector} <: BraidedSector
    a::I
end

sector_rank(::Type{TimeReversed{I}}) where {I<:BraidedSector} = sector_rank(I)

FusionStyle(::Type{TimeReversed{I}}) where {I<:BraidedSector}  = FusionStyle(I)
BraidingStyle(::Type{TimeReversed{I}}) where {I<:BraidedSector}  = BraidingStyle(I)
is_modular(::Type{TimeReversed{I}}) where {I<:BraidedSector}  = is_modular(I)
Nsymbol(a::TimeReversed{I}, b::TimeReversed{I}, c::TimeReversed{I}) where {I<:BraidedSector} = Nsymbol(a.a, b.a, c.a)
Fsymbol(a::TimeReversed{I}, b::TimeReversed{I}, c::TimeReversed{I}, d::TimeReversed{I}, e::TimeReversed{I}, f::TimeReversed{I}) where {I<:BraidedSector} = Fsymbol(a.a, b.a, c.a, d.a, e.a, f.a)
Rsymbol(a::TimeReversed{I}, b::TimeReversed{I}, c::TimeReversed{I}) where {I<:BraidedSector} = Nsymbol(a.a, b.a, c.a) * inv(Rsymbol(a.a, b.a, c.a))

Base.one(::Type{TimeReversed{I}}) where {I<:BraidedSector} = TimeReversed{I}(one(I))
Base.conj(c::TimeReversed{I}) where {I<:BraidedSector} = TimeReversed{I}(conj(c.a))
⊗(c1::TimeReversed{I}, c2::TimeReversed{I}) where {I<:BraidedSector} = map(TimeReversed{I}, c1.a ⊗ c2.a)
Base.IteratorSize(::Type{SectorValues{TimeReversed{I}}}) where {I<:BraidedSector} = Base.IteratorSize(SectorValues{I}())
Base.length(::SectorValues{TimeReversed{I}}) where {I<:BraidedSector} = Base.length(SectorValues{I}())
Base.getindex(::SectorValues{TimeReversed{I}}, i::Int) where {I<:BraidedSector} = TimeReversed{I}(SectorValues{I}()[i])
function Base.iterate(::SectorValues{TimeReversed{I}}, i::Int=0)  where {I<:BraidedSector}
    obj, next = iterate(SectorValues{I}(), i)
    return TimeReversed{I}(obj), next
end
findindex(::SectorValues{TimeReversed{I}}, a::TimeReversed{I})  where {I<:BraidedSector} = findindex(SectorValues{I}(), a.a)

Base.isless(c1::TimeReversed{I}, c2::TimeReversed{I}) where {I<:BraidedSector} = isless(c1.a, c2.a)


function S_matrix(::Type{𝒞}) where {𝒞<:BraidedSector}
    N = sector_rank(𝒞)
    S = zeros(ComplexF64, N, N)

    for i in 1:N, j in 1:N
        a = SectorValues{𝒞}()[i]
        b = SectorValues{𝒞}()[j]
        S[i, j] = sum(tr(Rsymbol(b, a, c) * Rsymbol(a, b, c)) * dim(c) for c in a ⊗ b)
    end

    return S
end

function T_vector(::Type{𝒞}) where {𝒞<:BraidedSector}
    N = sector_rank(𝒞)
    T = zeros(ComplexF64, N)

    for i in 1:N
        a = SectorValues{𝒞}()[i]
        T[i] = twist(a)
    end

    return T
end