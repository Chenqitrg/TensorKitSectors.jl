Z2grading(a::Irr{𝒞}) where {𝒞<:TY} = a == Irr{𝒞}(:σ) ? 1 : 0
rank(::Type{𝒞}) where {𝒞<:TY} = order(𝒞.parameters[1]) + 1

FusionStyle(::Type{Irr{𝒞}}) where {𝒞<:TY}  = SimpleFusion()
BraidingStyle(::Type{Irr{𝒞}}) where {𝒞<:TY}  = NoBraiding()
function Nsymbol(a::Irr{𝒞}, b::Irr{𝒞}, c::Irr{𝒞}) where {𝒞<:TY}
    if Z2grading(a)==Z2grading(b)==0
        return a.value * b.value == c.value
    elseif (Z2grading(a)==1 && Z2grading(b)==0)||(Z2grading(b)==1 && Z2grading(a)==0)
        return Z2grading(c)==1
    elseif Z2grading(a)==1 && Z2grading(b)==1
        return Z2grading(c)==0
    else
        throw(ArgumentError("Illegal simple objects $a, $b, $c"))
    end
end

function Fsymbol(a::Irr{𝒞}, b::Irr{𝒞}, c::Irr{𝒞}, d::Irr{𝒞}, e::Irr{𝒞}, f::Irr{𝒞}) where {𝒞<:TY}
    A, χ, ϵ = 𝒞.parameters
    is_match = Nsymbol(a, b, e) * Nsymbol(e, c, d) * Nsymbol(b, c, f) * Nsymbol(a, f, d)
    if Z2grading(a)==Z2grading(c)==0&&Z2grading(b)==1
        return χ(a.value, b.value)*is_match
    elseif Z2grading(a)==Z2grading(c)==1&&Z2grading(b)==0
        return χ(d.value, b.value)*is_match
    elseif Z2grading(a)==Z2grading(b)==Z2grading(c)==1
        return ϵ / sqrt(order(A)) / χ(e.value, f.value) * is_match
    else
        return is_match
    end
end

function Base.one(::Type{Irr{𝒞}}) where {𝒞<:TY}
    A = 𝒞.parameters[1]
    return Irr{𝒞}(identity_element(A))
end

Base.conj(c::Irr{𝒞}) where {𝒞<:TY} = Z2grading(c)==0 ? Irr{𝒞}(inverse(c.value)) : c
function ⊗(c1::Irr{𝒞}, c2::Irr{𝒞}) where {𝒞<:TY}
    if Z2grading(c1) == Z2grading(c2) == 0
        return (Irr{𝒞}(c1.value * c2.value),)
    elseif (Z2grading(c1) == 0 && Z2grading(c2) == 1)||(Z2grading(c1) == 1 && Z2grading(c2) == 0)
        return (Irr{𝒞}(:σ),)
    elseif Z2grading(c1) == Z2grading(c2) == 1
        return elements(𝒞.parameters[1])
    end
end

Base.IteratorSize(::Type{SectorValues{Irr{𝒞}}}) where {𝒞<:TY} = HasLength()
Base.length(::SectorValues{Irr{𝒞}}) where {𝒞<:TY} = rank(𝒞)

function Base.getindex(::SectorValues{Irr{𝒞}}, i::Int) where {𝒞<:TY}
    A = 𝒞.parameters[1]
    if i in 1:order(A)
        return Irr{𝒞}(𝒞.parameters[1][i])
    else
        return Irr{𝒞}(:σ)
    end
end

Base.iterate(::SectorValues{Irr{𝒞}}, i::Int=0)  where {𝒞<:TY} = i == rank(𝒞) ? nothing : (Irr{𝒞}[i], i + 1)

findindex(::SectorValues{Irr{𝒞}}, a::Irr{𝒞})  where {𝒞<:TY} = Z2grading(a)==0 ? findindex(a.value) : rank(𝒞)
Base.isless(c1::Irr{𝒞}, c2::Irr{𝒞}) where {𝒞<:TY} = isless(findindex(𝒞, c1), findindex(𝒞, c2))