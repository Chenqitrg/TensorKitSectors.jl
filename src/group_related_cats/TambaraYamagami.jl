struct TYIrr{A<:Group, χ, ϵ} <: FusionSector
    value::Any
    function TYIrr{A, χ, ϵ}(obj::Any) where{A<:Group, χ, ϵ}
        if !is_abelian(A)
            throw(ArgumentError("The group $A must be Abelian"))
        end
        if isa(obj, A) || obj == :σ
            new{A,χ,ϵ}(obj)
        else
            throw(ArgumentError("Illegal object $obj"))
        end
    end
end

Z2grading(a::TYIrr{A,χ,ϵ}) where {A<:Group,χ,ϵ} = a == TYIrr{A,χ,ϵ}(:σ) ? 1 : 0
sector_rank(::Type{TYIrr{A,χ,ϵ}}) where {A<:Group,χ,ϵ} = order(A) + 1

FusionStyle(::Type{TYIrr{A,χ,ϵ}}) where {A<:Group,χ,ϵ}  = SimpleFusion()
BraidingStyle(::Type{TYIrr{A,χ,ϵ}}) where {A<:Group,χ,ϵ}  = NoBraiding()
is_modular(::Type{TYIrr{A,χ,ϵ}}) where {A<:Group,χ,ϵ}  = false

function Nsymbol(a::TYIrr{A,χ,ϵ}, b::TYIrr{A,χ,ϵ}, c::TYIrr{A,χ,ϵ}) where {A<:Group,χ,ϵ}
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

function Fsymbol(a::TYIrr{A,χ,ϵ}, b::TYIrr{A,χ,ϵ}, c::TYIrr{A,χ,ϵ}, d::TYIrr{A,χ,ϵ}, e::TYIrr{A,χ,ϵ}, f::TYIrr{A,χ,ϵ}) where {A<:Group,χ,ϵ}
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

function Base.one(::Type{TYIrr{A,χ,ϵ}}) where {A<:Group,χ,ϵ}
    return TYIrr{A,χ,ϵ}(identity_element(A))
end

Base.conj(c::TYIrr{A,χ,ϵ}) where {A<:Group,χ,ϵ} = Z2grading(c)==0 ? TYIrr{A,χ,ϵ}(inverse(c.value)) : c
function ⊗(c1::TYIrr{A,χ,ϵ}, c2::TYIrr{A,χ,ϵ}) where {A<:Group,χ,ϵ}
    if Z2grading(c1) == Z2grading(c2) == 0
        return (TYIrr{A,χ,ϵ}(c1.value * c2.value),)
    elseif (Z2grading(c1) == 0 && Z2grading(c2) == 1)||(Z2grading(c1) == 1 && Z2grading(c2) == 0)
        return (TYIrr{A,χ,ϵ}(:σ),)
    elseif Z2grading(c1) == Z2grading(c2) == 1
        groupeles = elements(TYIrr{A,χ,ϵ}.parameters[1])
        return map(x->TYIrr{A,χ,ϵ}(x), groupeles)
    end
end

Base.IteratorSize(::Type{SectorValues{TYIrr{A,χ,ϵ}}}) where {A<:Group,χ,ϵ} = HasLength()
Base.length(::SectorValues{TYIrr{A,χ,ϵ}}) where {A<:Group,χ,ϵ} = rank(TYIrr{A,χ,ϵ})

function Base.getindex(::SectorValues{TYIrr{A,χ,ϵ}}, i::Int) where {A<:Group,χ,ϵ}
    if i in 1:order(A)
        return TYIrr{A,χ,ϵ}(TYIrr{A,χ,ϵ}.parameters[1][i])
    else
        return TYIrr{A,χ,ϵ}(:σ)
    end
end

Base.iterate(::SectorValues{TYIrr{A,χ,ϵ}}, i::Int=1)  where {A<:Group,χ,ϵ} = i > rank(TYIrr{A,χ,ϵ}) ? nothing : (SectorValues{TYIrr{A,χ,ϵ}}()[i], i + 1)

findindex(::SectorValues{TYIrr{A,χ,ϵ}}, a::TYIrr{A,χ,ϵ})  where {A<:Group,χ,ϵ} = Z2grading(a)==0 ? findindex(a.value) : rank(TYIrr{A,χ,ϵ})
Base.isless(c1::TYIrr{A,χ,ϵ}, c2::TYIrr{A,χ,ϵ}) where {A<:Group,χ,ϵ} = isless(findindex(SectorValues{TYIrr{A,χ,ϵ}}(), c1), findindex(SectorValues{TYIrr{A,χ,ϵ}}(), c2))