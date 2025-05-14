using TensorKit
using TensorKitSectors
abstract type BraidedCategory <: FusionCategory end
abstract type ModularCategory <: BraidedCategory end
abstract type ℨ{𝒞} <: ModularCategory end

function Irr{ℨ{𝒞}}(obj::GradedSpace{Irr{𝒞}, NTuple{N, Int}}, halfbraiding) where {𝒞<:FusionCategory, N}
    return Irr{ℨ{𝒞}}((obj, halfbraiding))
end

A = ℤ₂×ℤ₂
function χ(u::GroupElement{A}, v::GroupElement{A})
    u1, u2 = u.value
    v1, v2 = v.value
    sign = (-1)^(u1.value * v2.value + u2.value * v1.value)
    return sign
end
ϵ = 1
𝒟 = TY{A, χ, ϵ}
V = Vect[Irr{𝒟}](Irr{𝒟}(GroupElement{A}(GroupElement{ℤ₂}(0), GroupElement{ℤ₂}(0)))=>1, Irr{𝒟}(:σ)=>2)

T = rand(V⊗V←V⊗V)

Irr{ℨ{𝒟}}(V, T)