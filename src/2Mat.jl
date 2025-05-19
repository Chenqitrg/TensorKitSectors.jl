struct 𝔐at{N} <: FusionSector
    r::Int
    c::Int
    function 𝔐at{N}(r::Int, c::Int) where {N}
        if N<1
            throw(ArgumentError("$N should be larger than 1"))
        end
        if r<1 || r>N || c<1 || c>N
            throw(ArgumentError("The range of 1<$r<$N and 1<$c<$N condition is not satisfied"))
        end
        new{N}(r, c)
    end
end

sector_rank(::Type{𝔐at{N}}) where {N} = N^2

FusionStyle(::Type{𝔐at{N}}) where {N}  = SimpleFusion()
BraidingStyle(::Type{𝔐at{N}}) where {N} = NoBraiding()
Nsymbol(a::𝔐at{N}, b::𝔐at{N}, c::𝔐at{N}) where {N} = (a.c == b.r)&&(a.r==c.r)&&(b.c==c.c)
Fsymbol(a::𝔐at{N}, b::𝔐at{N}, c::𝔐at{N}, d::𝔐at{N}, e::𝔐at{N}, f::𝔐at{N}) where {N} = Nsymbol(a, b, e) * Nsymbol(e, c, d) * Nsymbol(b, c, f) * Nsymbol(a, f, d)

Base.conj(c::𝔐at{N}) where {N} = 𝔐at{N}(c.c, c.r)

⊗(c1::𝔐at{N}, c2::𝔐at{N}) where {N} = c1.c==c2.r ? (𝔐at{N}(c1.r,c2.c),) : ()
Base.IteratorSize(::Type{SectorValues{𝔐at{N}}}) where {N} = HasLength()
Base.length(::SectorValues{𝔐at{N}}) where {N} = N^2
Base.getindex(::SectorValues{𝔐at{N}}, i::Int) where {N} = 𝔐at{N}((i-1)÷N+1, (i-1)%N+1)
Base.iterate(::SectorValues{𝔐at{N}}, i::Int=1)  where {N} = i == N^2 ? nothing : (𝔐at{N}((i-1)÷N+1, (i-1)%N+1), i + 1)
findindex(::SectorValues{𝔐at{N}}, m::𝔐at{N})  where {N} = (m.r-1)*N + m.c
Base.isless(c1::𝔐at{N}, c2::𝔐at{N}) where {N} = isless(findindex(c1), findindex(c2))